# Parallélisation Level 2 de la reconstruction (surfacique + volumique)

> Reconstruction du maillage clippé avec **n procs pour le maillage A, m procs
> pour le maillage B, 1 proc clipper**, chaque maillage vivant dans son propre
> sous-communicateur (`FSClac` issu de `DivideIntoGroups`).

---

## Introduction — limites du code séquentiel

Le code séquentiel supposait implicitement **un seul proc par maillage**. Cette
hypothèse se cachait dans deux endroits qui ne se voyaient pas avec 1 proc mais
cassaient dès qu'un maillage était réparti sur plusieurs procs.

### 1. Numérotation locale des nœuds → collisions et trous

`FSCell2NodeBuilder` numérote les nœuds `0..N-1` **localement sur chaque proc**,
et `FSMeshReconstruction` écrivait cette numérotation locale comme
`GlobalNumber` (via `InitGlobalCellNumber(CT_Node, 0)` — offset 0 partout).

Or pour FSDM, `GlobalNumber` est le **seul identifiant stable inter-procs** :
au moment de recoller un maillage distribué, deux nœuds de même `GlobalNumber`
sont *le même nœud*. Conséquences avec ≥ 2 procs par maillage :

- **Collisions** : le nœud 5 du proc 1 et le nœud 5 du proc 2 sont des points
  physiquement différents mais portent le même numéro → FSDM les fusionne →
  mailles dégénérées.
- **Nœuds d'interface dupliqués** : un point sur la frontière entre deux
  partitions reçoit un numéro différent sur chaque proc → jamais recollé →
  trous dans le maillage.

Ce qui rendait le bug invisible en séquentiel : avec 1 proc, indices locaux ==
numéros globaux par coïncidence.

**Pourquoi la numérotation d'origine ne suffit pas** : la surface clippée
mélange trois populations de points — coins de A, coins de B, et **points
d'intersection nouvellement créés par le clipping** (Sutherland–Hodgman) qui
n'existent dans aucun des deux maillages d'entrée. Il faut donc une numérotation
neuve qui couvre les trois et soit cohérente entre procs. Le seul acteur qui
voit tous ces points simultanément est le **clipper**.

### 2. Prédicats géométriques non homogènes en dimension

Trois prédicats comparaient des quantités de dimension **longueur²** (produits
vectoriels, aires) à la tolérance qui est une **longueur**. Sur un maillage
unitaire c'est invisible ; sur un maillage à petite échelle (mailles ~2·10⁻⁴),
`cross ≈ 2·10⁻⁴ × 5·10⁻⁵ = 10⁻⁸ = tol` : tout se confond, des intersections
valides sont rejetées, des faces déclarées non couvertes à tort.

### 3. Cellules ghost traitées comme des cellules réelles

En numérotation locale FSDM, un proc détient ses cellules **owned** plus une
couronne de cellules **ghost** (répliques des cellules voisines). Le code
séquentiel émettait toutes les cellules sans distinction. Avec plusieurs procs :
faces dupliquées d'un proc à l'autre, et nœuds ghost référencés mais jamais
enregistrés.

---

## Architecture de la solution

Le clipper (proc 0, seul dans son groupe `MPI_COMM_SELF`) reçoit les faces des
deux maillages, calcule les matches, **assigne une numérotation globale**, et
renvoie à chaque proc ses matches enrichis. Chaque proc maillage reconstruit sa
part et **assemble la numérotation distribuée** conforme au contrat FSDM.

```
 procs maillage A          clipper (proc 0)          procs maillage B
 ────────────────          ────────────────          ────────────────
 extract faces ──GatherSend──▶ GatherReceiveAll
                              ComputeMatches
                              AssignGlobalNodeIds  (A)
      ◀──ScatterSend──────── (matches + IDs)
                              ComputeInvertedMatches
                              AssignGlobalNodeIds  (B)
                                        ─────ScatterSend────▶ (matches + IDs)
 BuildTopology                                          BuildTopology
 Reconstruct                                            Reconstruct
   └─ RemapToDistributedNumbering (AllGather dans le sous-comm du maillage)
```

---

## Étape 1 — Reconstruction surfacique

*(commits `56bd067` puis `791439a`)*

### 1a. Le clipper assigne la numérotation globale (côté production)

| Fichier | Modification |
|---------|-------------|
| `include/FSClipping/FSFaceMatcher.h` | `FSFaceMatch` gagne deux champs : `nodeGlobalIds` (un ID par point de `clippedPoly3D`) et `globalCellId` (un ID par match = futur ID de cellule Poly2D). Déclaration de `AssignGlobalNodeIds`. |
| `src/FSFaceMatcher.cpp` | **`AssignGlobalNodeIds(matches, tol)`** : déduplication géométrique par `NodeKey` (quantification par tolérance), IDs contigus `0..K-1`. Deux matches partageant un point géométrique reçoivent le **même** ID, même s'ils partent vers des procs différents. `GetBufSize`/`Pack`/`Unpack` étendus pour sérialiser les nouveaux champs. |
| `src/FSClippingEngine.cpp` | `RunMatcherProc` appelle `AssignGlobalNodeIds` pour chaque côté (A, puis B après `ComputeInvertedMatches`) juste avant le `ScatterSend`. |

### 1b. La reconstruction consomme la numérotation (côté consommation)

| Fichier | Modification |
|---------|-------------|
| `include/FSClipping/FSCell2NodeBuilder.h` | Surcharge `AddClippedPolygon(cellId, poly, globalIds)` + `NodeGlobalNumbers()` + map interne `keyToGlobalId_`. |
| `src/FSCell2NodeBuilder.cpp` | La surcharge mémorise `NodeKey → ID clipper`. `NodeGlobalNumbers()` retourne, pour chaque nœud local, son ID clipper (ou -1). Correction annexe d'un vrai bug : `"...cellId : " + globalId` faisait de l'arithmétique de pointeur au lieu d'une concaténation. |
| `include/FSClipping/FSTopologyAssembler.h` | `FSTopologyData` gagne `nodeGlobalNumbers` (parallèle à `globalCoords`) et `poly2DGlobalNumbers` (parallèle aux cellules Poly2D). |
| `src/FSTopologyAssembler.cpp` | `BuildSurfaceTopo` remplit ces champs **uniquement si** les matches portent des IDs (mode parallèle) ; sinon vides (mode séquentiel inchangé). |

### 1c. Prédicats géométriques scale-aware *(commit `791439a`)*

| Fichier | Correction |
|---------|-----------|
| `include/FSClipping/FSClippingUtil.h` | **`IsPointInsideEdge`** (le bug principal) : `cross >= -tol` → `cross >= -tol·\|B−A\|` (distance signée point–droite). Sans ça, les arêtes du clipper ne coupaient pas le sujet → sur-couverture → faces déclarées NOT_COVERED à tort → matches supprimés. **`ArePointsColinear2D`** : `\|cross\| > tol·\|dir\|`. Ajout de `PolygonPerimeter`. |
| `src/FSFaceMatcher.cpp` | **`ComputeMatch`** et **`PreserveUncoveredFaces`** : seuils d'aire exprimés en largeur moyenne `tol · périmètre/2` (dimension aire = longueur × longueur). |

**Effet mesuré** (test 2 cubes 4×4 / 5×5, 5 procs) : matches passés de 41/64 à
64/64, aire d'interface reconstruite de 2,31·10⁻⁶ (overlaps) à **exactement**
1,00·10⁻⁶.

### 1d. Assemblage distribué de la numérotation *(commit `791439a`)*

| Fichier | Modification |
|---------|-------------|
| `include/FSClipping/FSMeshReconstruction.h` | Déclaration de `RemapToDistributedNumbering`, ajout du membre `tol_`. |
| `src/FSMeshReconstruction.cpp` | Cœur de la solution — voir ci-dessous. |

**`RemapToDistributedNumbering`** applique le contrat FSDM (documenté dans
`FSMesh.h`) :

1. **AllGather** des listes de nœuds de chaque proc (clé géométrique = coord
   quantifiée par tol) sur le sous-communicateur du maillage.
2. **Ownership déterministe min-rank** : en scannant les procs par rang
   croissant, le premier qui liste un nœud le possède ; la numérotation en ordre
   de scan donne les **plages contiguës par proc** exigées par FSDM.
3. `InitUnstructNodes` reçoit **uniquement les nœuds owned** ; les nœuds
   d'interface n'existent plus qu'une fois.
4. **Remap in place** de la connectivité (`cell2NodePoly2D`, `cell2NodePoly3D`)
   vers les numéros distribués.
5. L'attribut `GlobalNumber` garde ces numéros comme identifiant stable.

Le mode 1-proc-par-maillage passe par un chemin identité, inchangé.

**Point subtil** : `polyFaces` (stockage des faces Poly3D) utilise des indices
locaux à la cellule et **ne doit PAS** être remappé.

---

## Étape 2 — Reconstruction volumique

*(modifications actuelles, non commitées au moment de l'écriture)*

Le volumique réutilise tout l'appareillage de l'étape 1 (numérotation clipper,
`RemapToDistributedNumbering` qui remappe déjà `cell2NodeInner`). Deux
problèmes restaient, tous deux liés aux **cellules ghost**.

| Fichier | Modification |
|---------|-------------|
| `src/FSTopologyAssembler.cpp` | **`BuildVolumeTopo`** : skip des cellules ghost (`if(!cellPool.IsOwned(c)) continue`) lors de l'ajout des cellules volumiques — sinon les cellules non-owned, répliquées chez les voisins, seraient émises deux fois. La boucle de `cellParent` itère désormais sur le nombre réel de cellules du builder (owned) et non sur le compte du pool. |
| `src/FSTopologyAssembler.cpp` | **`UpdateOldCell2Node`** (surfaces marker non clippées) reçoit maintenant le `FSCellPool` et **ne traite que les cellules owned**. La taille du tableau de sortie est calculée sur le compte réel des cellules gardées (owned ∧ non clippées) et non sur `nCellOld − bdry.size()`. C'est le pendant surfacique du skip ghost de `BuildVolumeTopo`. |
| `src/FSTopologyAssembler.cpp` | **`AppendUnclippedSurfaces`** passe le `cellPool` à `UpdateOldCell2Node` ; include de `<FSMeshData.h>`. |
| `src/FSClippingEngine.cpp` | `BuildTopology` mode Volume : pool vide (`kEmptyPool`) pour les types de cellules sans face d'interface (proc dont la partition ne touche pas le marker) ; propagation de `poly2DGlobalNumbers` du surfacique au volumique. |
| `src/FSMeshReconstruction.cpp` | Le remap est activé dès `nProcs > 1` (fusion géométrique générique, valable nœuds de surface **et** nœuds volumiques). `cell2NodeInner` est remappé au même titre que les registres Poly2D/Poly3D. |

### Le bug de deadlock résolu

Symptôme : le test volumique 5 procs faisait un **deadlock MPI**.

Cause : `UpdateOldCell2Node` traitait toutes les cellules ghost incluses. Comme
`BuildVolumeTopo` skippait déjà les ghosts, les nœuds exclusifs aux cellules
ghost n'étaient jamais enregistrés dans le builder → `coord2Node.at(key)`
**lançait une exception** sur les procs de rang local 0 de chaque groupe.
L'exception non catchée les bloquait **avant** l'`AllGather` collectif de
`RemapToDistributedNumbering`, pendant que les autres procs du groupe
attendaient indéfiniment sur ce collectif.

Fix : le skip des ghosts dans `UpdateOldCell2Node` élimine les nœuds manquants,
donc l'exception, donc le deadlock.

---

## Étape 2bis — Robustesse aux partitions fines (test 2 cylindres)

Le passage des cubes (grille structurée) aux **2 cylindres** — maillages plus
gros, partitionnés plus finement — a exposé deux bugs qui ne se déclenchaient
pas sur les cubes, car ils dépendent de la topologie des partitions produites
par RCB. Testé jusqu'à **9 procs** (4 procs par maillage + 1 clipper).

| Fichier | Modification |
|---------|-------------|
| `include/FSClipping/FSPolyFaceBuilder.h`, `src/FSPolyFaceBuilder.cpp` | Ajout de **`HasCell(globalId)`** : indique si une cellule (par son id global) est connue du builder de faces, c.-à-d. porte un match de surface clippée. |
| `src/FSTopologyAssembler.cpp` | **`BuildVolumeTopo`** : la boucle `AddInnerFaces` ne traite que les cellules `IsOwned(c) && polyFaceBuilder_->HasCell(c)`. |
| `src/FSTopologyAssembler.cpp` | **`BuildSurfaceTopo`** : le cas `matches` vide ne fait plus `SetAndPrintAndExit`, il construit une topologie de surface vide mais valide. |

### Bug 1 — cellules frontière sans match

`FSCell2NodeBuilder : Unknown cellId` dans `AddInnerFaces` → `LocalCellIndex`.

`bdryCellPool` (les cellules ayant une face sur le marker, issues de
l'extraction géométrique **locale**) peut contenir des cellules **sans match** :
leur face marker n'a pas de correspondance dans l'autre maillage, ou le match a
été routé vers un autre proc. Ces cellules ne sont pas dans le
`polyFaceBuilder_` (qui ne connaît que les cellules matchées, `elemOwner1`), donc
`AddInnerFaces → LocalCellIndex` échouait.

Trace confirmant le diagnostic : `Unknown cellId=4386, cellId2L_.size=11` — le
builder ne connaît que 11 cellules matchées mais le pool en réclame une
douzième.

**Pourquoi c'est cohérent de la skipper** : une cellule frontière sans match
n'a pas reçu de polygone clippé → elle reste `onTheBorder = false` → elle est
classée comme cellule volumique normale (Hexa), pas Poly3D. Elle n'a donc aucune
face interne à reconstruire. Le skip via `HasCell` s'aligne exactement sur ce
classement.

### Bug 2 — robustesse au proc sans face d'interface

`BuildSurfaceTopo` faisait `SetAndPrintAndExit("No matches were found")` dès
qu'un proc n'avait aucun match. Un sous-domaine ne touchant pas la frontière
(0 face extraite → 0 match) **crasherait** donc, provoquant un deadlock des
autres procs du groupe sur le collectif `AllGather` de reconstruction.

Ce cas n'est pas apparu dans les tests (RCB répartit la frontière du cylindre
sur **tous** les sous-domaines — comptage vérifié : 133–143 faces/proc à 9
procs, aucun proc à 0), mais le code est désormais **robuste** : un proc sans
match construit une topologie de surface vide et participe normalement à la
reconstruction collective ; l'étape volumique passe toutes ses cellules comme
volume non clippé.

En mode **séquentiel**, l'ancien message d'erreur « No matches were found » (qui
signalait un mauvais marker) disparaît : un clipping sans match produit
simplement un maillage sans surface clippée. Aucun test ne dépendait de ce
message.

---

## Tests

| Fichier | Contenu |
|---------|---------|
| `test/FSClipping/FSClippingTestGlobalNumbering.cpp` *(commit `56bd067`)* | Numérotation clipper : déduplication, tolérance, contiguïté, propagation `AssignGlobalNodeIds` → `BuildSurfaceTopo`, non-régression du mode séquentiel. |
| `test/FSClipping/Parallel/FSClippingTestInterfaceParallel.cpp` | Tests `SurfaceInterface` et `VolumeInterface` (cubes, 5 procs), et **`VolumeInterface2Cylinders`** (cylindres, 5 et 9 procs), tous avec export partition-independent qui valide le recollement distribué via `GlobalNumber`. |

### Validation finale

| Test | Config | Résultat |
|------|--------|----------|
| `SurfaceInterface` (cubes) | 5 procs | PASSED |
| `VolumeInterface` (cubes) | 5 procs | PASSED |
| `VolumeInterface2Cylinders` | 5 procs | PASSED |
| `VolumeInterface2Cylinders` | 9 procs | PASSED |
| Suite séquentielle | 1 proc | 35/35 |

- **Export partition-independent** (recollement distribué via `GlobalNumber`) :
  réussi — preuve que le contrat FSDM est respecté.
- **Analyse des maillages recollés (cubes)** : 0 incohérence
  `GlobalNumber`↔coordonnées dans les deux sens, 0 hexa dégénéré, 0 cellule
  dupliquée. Volume hexa conservé exactement (coarse 48/64 = 7,5·10⁻¹⁰ ;
  fine 100/125 = 8·10⁻¹⁰).
- **Extraction équilibrée (cylindres)** : chaque proc extrait des faces
  d'interface, aucun sous-domaine vide (266–285 faces/proc à 5 procs, 133–143 à
  9 procs, 77–80 à 15 procs).
- **Non-régression** : 35/35 tests séquentiels ; `FSClippingTestMatchesPar`
  inchangé à 3 procs (`Identical` PASSED, `Intersection`/`Included` en échec
  **préexistant**).

---

## Étape 2ter — Répartition asymétrique (n ≠ m procs par maillage)

Jusqu'ici les deux maillages recevaient le **même** nombre de procs (parité de
rang : pairs → maillage 1, impairs → maillage 2). Test d'une répartition
**asymétrique** : 4 procs pour le cylindre 1, 2 procs pour le cylindre 2, 1
clipper (7 procs), via un `meshID` par **plages de rangs** contiguës
(`AsymMeshID`) au lieu de la parité.

Le routage clipper (`FSFaceExchange::GatherReceiveAll` / `ScatterSend`) est déjà
**agnostique au layout** : il reconstruit `localToGlobal` / `cellToGlobalProc` à
partir des en-têtes (`worldProcID`, `localProcID`, `meshID`), sans hypothèse de
parité. Vérifié : le clipper voit bien 4 procs côté A et 2 côté B, 551+551 faces,
2906 matches. **L'asymétrie ne casse donc pas le clipping.**

### Bug exposé (corrigé)

L'asymétrie rend fréquent le cas d'un proc à **partition clippée vide** (après
attribution min-rank des nœuds, le localRank 0 du cylindre 1 se retrouve avec 0
nœud possédé et 0 Poly2D/Poly3D).

| Fichier | Correction |
|---|---|
| `src/FSMeshReconstruction.cpp` | **`CopyCellAttributes`** : `if(numCells == 0) continue;` en tête de boucle. Sans ce garde, `topo.cellParentType.at(polyType)[0]` déréférençait un tableau vide → **segfault** (backtrace `CopyCellAttributes` ← `Build` ← `Reconstruct`). Un `FSIntArrayT` de taille 0 passé à `InitCellAttribute` déclenche aussi son contrôle `NDims()`. |

Précisions importantes vérifiées dans les sources FSDM :
- `InitCellAttribute(name, values)` est **purement local** (adopte le tableau,
  contrôle de taille vs `GetNCells`, aucun collectif MPI) — donc sauter un pool
  vide sur certains procs est sûr. Le « deadlock collectif » soupçonné
  précédemment était une fausse piste ; le vrai déclencheur était ce
  segfault/tableau-vide.
- Un proc possédant 0 nœud est **légal** : `InitNodePool(nNodes)` fait un exscan
  et n'abort que si `node2Proc[nProcs] == 0`, c.-à-d. si **tout** le
  sous-communicateur a 0 nœud.
- Piège de test : le **marqueur de bord est un attribut** (`CADGroupID`, lu dans
  `FSFaceSeparator`). Un « strip all attributes » naïf le supprime → `Extract`
  renvoie 0 face → tout le groupe reconstruit du vide. `StripCellAttributes`
  conserve `CADGroupID` + `GlobalNumber` (`AT_CADGroupID` est dans `FSEnums`).

### Validation

- `SurfaceInterface2CylindersAsymNoAttr` (7 procs, 4+2+1) : **PASSED**.
- `VolumeInterface2CylindersAsymNoAttr` (7 procs, 4+2+1) : **PASSED** (plus de
  segfault).
- Non-régression : séquentiel 35/35 ; symétrique 5 procs (teardown export/import
  complet) exit 0.

### Reste ouvert (hors clipping)

L'**export/import HDF5 partition-independent** du teardown (`polyMeshExportImport`)
**hang** quand un proc du groupe détient une partition poly vide. C'est une
divergence collective dans la couche I/O de FSDM, en aval de la reconstruction
(qui est correcte). Les tests asym restreignent leur teardown à `CheckMesh` en
attendant. À traiter soit en garantissant qu'aucun proc ne reçoit une partition
clippée vide, soit côté FSDM.

---

## Lancer les tests parallèles (3 / 5 / 7 procs)

Chaque test parallèle est écrit pour un **layout de procs précis** et se
**skippe** proprement (`GTEST_SKIP`, via `SKIP_UNLESS_EXACT_PROCS` /
`SKIP_UNLESS_MIN_PROCS` dans `TestUtilsParallel.hpp`) quand le nombre de procs
courant ne correspond pas — plus aucun crash/deadlock quand on lance au mauvais
nombre.

| Test | Layout | Procs (CTest) |
|---|---|---|
| `FSClippingTestMatchesPar.*` | 3 fixes | exactement 3 |
| `SurfaceInterface` / `VolumeInterface` (cubes) | parité | exactement 3 |
| `VolumeInterface2Cylinders` (symétrique) | parité | exactement 5 — 9 en manuel |
| `SurfaceInterface2CylindersAsymNoAttr` / `VolumeInterface2CylindersAsymNoAttr` | 4+2+1 | exactement 7 |

Les tests pinés à un proc-count exact le sont parce que leur teardown
export/import HDF5 déclenche le hang I/O à haut nombre de procs (> ~5, cf.
« Limitations ») ; on les exécute donc uniquement là où ils terminent. Les tests
asym (7 procs) n'ont pas ce problème : leur teardown est restreint à `CheckMesh`
(pas d'export/import).

Le même code de test est enregistré dans CTest à **3, 5 et 7 procs** (une
`OBJECT` library partagée, un exécutable par proc-count :
`FSClippingParallelTest{3,5,7}`, préfixes CTest `np3.` / `np5.` / `np7.`). Ainsi
**un simple « Run all tests » dans VSCode** exécute les trois configurations, et
chaque test tourne uniquement là où il a un sens (les autres `SKIPPED`).

Lancement manuel d'un proc-count donné :

```bash
mpirun -np 3 ./test/FSClippingParallelTest3   # match send/receive + interface cubes
mpirun -np 5 ./test/FSClippingParallelTest5   # + cylindres symétrique
mpirun -np 7 ./test/FSClippingParallelTest7   # + cylindres asymétrique 4+2+1
```

---

## Limitations connues / suite

- Les tests parallèles `FSClippingTestMatchesPar.Intersection` et `.Included`
  (3 procs) échouent — **préexistant**, indépendant de ce chantier.
- Le cas où plusieurs types d'éléments volumiques touchent la frontière : dans
  `BuildVolumeTopo`, le dernier type écrase le `cellParent[CT_Poly3D]`
  (voir `NOTE` dans le code).
- Validé sur cubes hexaédriques structurés (5 procs) et 2 cylindres (5 et 9
  procs). Le cas d'un sous-domaine sans face d'interface est géré défensivement
  mais n'a pas encore été rencontré en pratique — à provoquer explicitement pour
  le couvrir par un test dédié. À étendre aux configurations rotor/stator.
- Test 15 procs : abort en *oversubscribe* (plus de procs que de cœurs
  physiques) — limite d'environnement MPI, pas du code ; l'extraction des faces
  se déroulait correctement sur les 14 procs maillage avant l'abort.
