# Parallelization Roadmap – FSClipping

## Contexte

L'algorithme de clipping nécessite les faces frontières des deux maillages.
En séquentiel, un seul processus possède tout. En parallèle, chaque proc
possède UNE partition du maillage. L'objectif est de distribuer le calcul
sans créer de goulot d'étranglement.

Configuration de référence :
- `globalClac` = MPI_COMM_WORLD (N procs au total)
- Chaque proc reçoit un `clac` à 1 proc via `DivideIntoGroups(meshID, clac)`
- Proc 0..P-1 : procs "mesh1"
- Proc P..Q-1 : procs "mesh2"
- Proc Q..N-1 : procs "clippers"

---

## Niveau 1 — Proc clipper centralisé

### Objectif

Un seul proc clipper (proc 2) reçoit toutes les faces des deux maillages
et effectue l'intégralité du calcul de clipping.

### Schéma

```
Proc 0 (mesh1) → PackFaces → MPI_Send → Proc 2
Proc 1 (mesh2) → PackFaces → MPI_Send → Proc 2
                                         Proc 2 : UnpackFaces
                                                  FSFaceMatcher
                                                  BuildSurfaceTopo
                                                  FSMeshReconstruction
```

### Sous-objectifs

- [x] **L1.1** `FSFaceSerializer::Pack / Unpack` — sérialisation géométrie-seule
       des `FSClippingFace` en buffer `double[]` pour MPI
- [x] **L1.1** Test unitaire séquentiel de Pack/Unpack (round-trip exact)
- [ ] **L1.2** Extraction par proc : chaque proc extrait ses faces frontières
       via `FSBoundaryFaceProvider::Extract` avec son `clac` local
- [ ] **L1.3** Envoi des faces au proc clipper via `MPI_Send / MPI_Recv`
       sur `globalClac`
- [ ] **L1.4** Calcul des matches sur le proc clipper avec `FSFaceMatcher`
       (clac local à 1 proc → BVH purement local)
- [ ] **L1.5** Test parallèle à 3 procs validant L1.2 → L1.4

### Limites

- Tout passe par un seul proc → goulot mémoire et calcul si N_faces >> 1
- Pas d'exploitation du parallélisme pour le matching lui-même

---

## Niveau 2 — K procs clippers avec partition spatiale

### Objectif

Distribuer le travail de matching sur K procs clippers en découpant l'espace
en K régions. Chaque clipper ne reçoit que les faces dont la boîte englobante
intersecte sa région.

### Schéma

```
Proc 0 (mesh1) : faces1_A, faces1_B       Proc 1 (mesh2) : faces2_A, faces2_B
        │                                          │
        ├── faces1_A ──────────────────────────────┤──→ Proc 2 (clipper A)
        └── faces1_B ──────────────────────────────┤──→ Proc 3 (clipper B)
                                                   │
                                           faces2_A──→ Proc 2
                                           faces2_B──→ Proc 3
```

### Sous-objectifs

- [ ] **L2.1** Calcul d'une partition spatiale globale (K boîtes englobantes)
       à partir des bounding boxes de toutes les faces — via `MPI_Allgather`
       des bounding boxes sur `globalClac`
- [ ] **L2.2** Routage : chaque proc mesh envoie chaque face au clipper
       dont la région intersecte la bounding box de la face
       (`FSBoundingBoxUtil` pour le test d'intersection)
- [ ] **L2.3** Chaque clipper reçoit ses faces (mesh1 + mesh2), exécute
       `FSFaceMatcher` localement
- [ ] **L2.4** Gestion des faces frontières de région (face intersectant
       deux régions → envoyée aux deux clippers concernés)
- [ ] **L2.5** Test parallèle à 4 procs (1 mesh1 + 1 mesh2 + 2 clippers)

### Limites

- Les faces à cheval sur plusieurs régions sont dupliquées → peut créer
  des doublons de matches qu'il faut dédupliquer
- Nécessite un BVH global des régions (coût de construction O(N log N))

---

## Niveau 3 — Pipeline totalement distribué N×M→K

### Objectif

Généraliser à N procs pour mesh1, M procs pour mesh2, K procs clippers.
Chaque proc mesh gère une partition locale et envoie ses faces aux clippers
concernés sans passer par un proc centralisé.

### Schéma

```
N procs mesh1 ──→ routage BVH-guidé ──→ K procs clippers
M procs mesh2 ──→ routage BVH-guidé ──→ K procs clippers
                                         résultats distribués sur K procs
```

### Sous-objectifs

- [ ] **L3.1** BVH global des régions clippers construit collectivement
       (chaque clipper annonce sa région, BVH diffusé à tous via `MPI_Bcast`)
- [ ] **L3.2** Chaque proc mesh1/mesh2 interroge le BVH global pour déterminer
       à quels clippers envoyer chacune de ses faces locales
       → utilisation de `MPI_Isend` non-bloquant
- [ ] **L3.3** Chaque clipper reçoit les faces en streaming
       (`MPI_Irecv` + `MPI_Waitany`) et construit son BVH local
       au fil des réceptions
- [ ] **L3.4** Renvoi des résultats (matches / maillage clippé)
       aux procs mesh1 d'origine pour reconstruction du maillage distribué
- [ ] **L3.5** Dépendance cyclique `FSMeshReconstruction` : adaptation pour
       recevoir une topologie distribuée sur K procs et reconstruire
       le maillage sur les N procs mesh1
- [ ] **L3.6** Tests parallèles à 6 procs (2+2+2)

### Limites

- Complexité d'implémentation élevée (communication non-bloquante, deadlock)
- Nécessite de rendre `FSMeshReconstruction` conscient de la distribution

---

## Récapitulatif

| Niveau | Procs | Matching distribué | Mémoire | Complexité |
|--------|-------|--------------------|---------|------------|
| L1     | 3     | Non (1 clipper)    | O(N)    | Faible     |
| L2     | 4+    | Partiel (K clippers, partition spatiale) | O(N/K) | Moyenne |
| L3     | 6+    | Total (N×M→K)      | O(N/K)  | Élevée     |
