# Extraction des faces de frontière en parallèle

> Correction de `FSFaceSeparator::SeparateFaces()` : le nombre de faces extraites
> à l'interface dépendait du nombre de procs. Alignement sur
> `FaceBasedMeshAdapter::SeparateFaces()` / `ExchangeSeparatedBdryFaces()` de CODA.

---

## Symptôme

Sur un cas rotor/stator de plusieurs dizaines de millions d'éléments, le nombre
total de faces extraites à l'interface changeait avec le découpage :

| maillage | 35 procs | 40 procs | écart |
|---|---|---|---|
| rotor  (`interface_rotor`)  | 159 211 | 158 996 | −215 |
| stator (`interface_stator`) | 152 737 | 152 413 | −324 |

Le clipping allait au bout, produisait deux maillages, et `mesh.check()` passait.
Les faces manquantes restaient de simples `Quad4` non clippés : la topologie
restait cohérente, seule la géométrie de l'interface était fausse. **Perte
silencieuse**, donc, et non un crash.

---

## Cause

### Comment FSDM apparie une face

`FSMeshFaceExtractor::PrepareFaceConnectivity()` insère les cellules dans un
ordre imposé — **pseudo → volume → surface** — et n'insère que les cellules
*owned*. La première cellule qui insère une face devient `mOwner`, la seconde
`mNeighbor`.

Pour une face d'interface, deux cellules se rencontrent : un `Hexa8` (volume) et
un `Quad4` (surface). Les volumes étant insérés en premier :

```
mOwner    = Hexa8   → mFaceIndex >= 0   (indice de la face dans l'hexaèdre)
mNeighbor = Quad4   → mFaceIndex == -1  (une cellule surface n'a pas d'indice de face)
```

Point essentiel : **seul le `Quad4` porte l'attribut `CADGroupID`**. L'`Hexa8`
ignore qu'il touche l'interface. Décider si une face doit être clippée impose
donc de lire l'attribut sur le `Quad4`.

### Le code avant correction

```cpp
if(face.mOwner.mFaceIndex >= 0) {              // owner est un volume
  if(face.mNeighbor.mCellProcID == procID) {   // <-- le bug
    if(IsUnstructSurfaceCellType(face.mNeighbor.mCellType)) {
      marker = boundaryMarkers(face.mNeighbor.mCell);   // lecture LOCALE
      if(marker == markerBoundaryClipped)
        bdryFaces.emplace_back(...);
    }
  }
}
```

`boundaryMarkers(face.mNeighbor.mCell)` indexe le pool de cellules **local**.
Le garde `mCellProcID == procID` est donc nécessaire : sans lui on lirait un
indice invalide. Le problème n'est pas le garde, mais l'absence de traitement
du cas qu'il écarte.

### Ce que fait le partitionneur

RCB partitionne les `Hexa8` et les `Quad4` **indépendamment**, par coordonnées.
Rien ne garantit qu'un `Quad4` reste sur le proc de l'`Hexa8` qu'il ferme.

Quand la paire est séparée entre le proc **P** (volume) et le proc **Q**
(surface), aucun appariement local n'est trouvé. `FindUmatchedFacesViaBroadcast()`
recolle alors les faces orphelines par diffusion — et la même face physique
apparaît **des deux côtés, owner et neighbor inversés** :

| | `mOwner` | `mNeighbor` |
|---|---|---|
| sur **P** (a l'`Hexa8`) | `Hexa8`, local | `Quad4`, sur Q |
| sur **Q** (a le `Quad4`) | `Quad4`, local | `Hexa8`, sur P |

Appliqué à l'ancien code :

- **sur P** : `mOwner.mFaceIndex >= 0` ✔, mais `mNeighbor.mCellProcID == Q ≠ P`
  → rejetée par le garde ;
- **sur Q** : `mOwner` est le `Quad4`, donc `mFaceIndex == -1`
  → rejetée par le tout premier `if`.

**Aucun proc ne conserve la face.**

Le nombre de paires ainsi séparées dépend entièrement du découpage, d'où la
variation avec le nombre de procs.

### Reproduction en miniature

Sur `cube_hexa_coarse_par.grid` (16 faces sur le marker 6), avec le code d'origine :

| procs | 1 | 2 | 4 | 5 | 7 |
|---|---|---|---|---|---|
| faces extraites | 16 | 16 | **15** | 16 | 16 |

### Pourquoi le bug est resté invisible si longtemps

Les tests sur `mesh_cylinder` à 5, 10 et 20 procs donnaient toujours 551 faces.
Mesures faites avec le code d'origine :

| procs | 5 | 10 | 20 | 40 |
|---|---|---|---|---|
| paires séparées par RCB | **0** | 4 | 15 | 42 |
| faces extraites | 551 | 551 | **550** | 551 |

Deux effets se combinent :

1. **À 5 procs, RCB ne sépare aucune paire** sur ce maillage — le bug ne peut
   pas se déclencher.
2. **Une paire séparée n'est perdue que si son `Quad4` porte le marker clippé.**
   `mesh_cylinder` compte 551 faces d'interface pour 1624 paires volume/surface
   (~1/3) : les paires séparées tombent le plus souvent sur un autre marker.
   D'où 42 paires séparées à 40 procs et pourtant le bon total.

Corollaire important : **un total correct sur un run ne prouve rien**. C'est ce
qui rendait ce bug difficile à isoler par simple comparaison de compteurs.

---

## Correction

L'approche est celle de CODA (`FaceBasedMeshAdapter`), qui traite ce cas depuis
toujours. Les deux moitiés d'une paire séparée y portent un nom :

| liste | vue depuis | possède | lui manque |
|---|---|---|---|
| `sepBdryFaces` (*separated*) | **Q**, proc du `Quad4` | le **marker**, lisible localement sur `mOwner` | la cellule volume, distante |
| `extBdryFaces` (*external*)  | **P**, proc de l'`Hexa8` | la **face du bon côté**, celle dont le clipping a besoin | le marker, illisible |

**Aucune des deux moitiés ne suffit seule** : l'une a l'information, l'autre a la
face utile. L'échange MPI les recolle.

Ces deux listes sont des **variables locales** à `SeparateFaces()` : elles
servent uniquement à porter l'échange et ne sortent jamais de la fonction. Dans
CODA elles sont retournées à l'appelant, car le solveur les réutilise (conditions
aux limites, liens de halo) ; ici le clipping n'a besoin que de `bdryFaces`.

### La boucle principale

```cpp
if(face.mOwner.mFaceIndex >= 0) {              // owner = volume
  if(!IsUnstructSurfaceCellType(face.mNeighbor.mCellType))
    continue;                                  // face intérieure
  if(face.mNeighbor.mCellProcID == procID)
    → bdryFaces      // paire locale complète : inchangé
  else
    → extBdryFaces   // le Quad4 est ailleurs : on garde la face, on attend le marker
} else {                                       // owner = surface
  → sepBdryFaces     // on lit le marker localement, on l'enverra
}
```

### L'échange

`ExchangeSeparatedBdryFaces()` suit le modèle **push** de CODA : le proc qui
possède le `Quad4` envoie spontanément le marker à celui qui possède l'`Hexa8`.

```
Q (a le marker)  --[ connectivité, marker, globalNumber ]-->  P (a la face utile)
                                                               └─> bdryFaces
```

La face conservée est celle de **P**, côté volume — exactement ce dont le
clipping a besoin. Seule l'information voyage.

Les asserts de vérification de CODA sont repris tels quels :

```cpp
assert(sbf.mOwner    == ebf._faceFSDM->mNeighbor);  // l'owner distant est mon voisin
assert(sbf.mNeighbor == ebf._faceFSDM->mOwner);     // et réciproquement
```

ainsi que le contrôle de cohérence global, qui vérifie que chaque `ext` a bien
son `sep` :

```cpp
auto missingBdryFaces = static_cast<FS_int32T>(extBdryFaces.size() - sepBdryFaces.size());
fsmesh.GetClac()->SumInPlace(&missingBdryFaces, 1);
assert(missingBdryFaces == 0);
```

### Absence de doublon

Une paire séparée est vue par deux procs, mais un seul la conserve :

- **P** (volume) la reçoit et la pousse dans `bdryFaces` ✔
- **Q** (surface) n'envoie que le marker, ne garde rien ✘

C'est ce que verrouille le test `NoFaceIsExtractedTwice`.

---

## Deux écarts par rapport à CODA

Ce ne sont pas des transpositions littérales, et ils comptent.

### 1. L'ordre de tri

CODA apparie les faces **par position** : les deux listes triées, la n-ième face
envoyée est la n-ième attendue.

Or l'`operator<` de FSDM ordonne par `(mOwner, mNeighbor)` — et owner/neighbor
sont **inversés** entre les deux moitiés. Trier chaque liste avec sa propre clé
ne donne donc pas le même ordre des deux côtés. Observé sur `np6.TwoCylindersVolume` :

```
proc 0 (EXT), trié par ownCell :  2967, 3151, 3335, 3394, ...
proc 2 (SEP), trié par ownCell :  2496, 2508, 2520, 2808, ...
```

L'appariement étant `ownCell(EXT) ↔ nbrCell(SEP)`, les deux suites ne sont pas
co-monotones (`2496→4014` casse l'ordre) : les faces s'appariaient de travers, et
les asserts de CODA — repris ici — l'ont immédiatement détecté.

Les deux moitiés sont donc triées sur **la même clé**, celle vue depuis le volume :

```cpp
(volume.mCellProcID, volume.mCell, surface.mCellProcID, surface.mCell)
```

### 2. Sends non bloquants

`FSClac::SendBuffer()` est **bloquant** : il attend que le destinataire poste son
`ReceiveBuffer()`. Une première version enchaînait tous les envois puis toutes les
réceptions — ce qui interbloque dès que deux procs détiennent chacun des faces
séparées pour l'autre : A bloque en envoi vers B pendant que B bloque en envoi
vers A, et aucun n'atteint sa boucle de réception.

Constaté à 35+35 procs : 70 procs entraient dans `Extract`, **33 en sortaient**.

CODA n'a pas ce problème car il emploie `PostSend`/`PostRecv`, **non bloquants**,
finalisés par `WaitCommReq`. Traduction retenue :

| CODA | FSClac |
|---|---|
| `PostSend` | `NBSendBuffer(proc, false)` |
| `PostRecv` | `NBReceiveBuffer(proc, false)` |
| `WaitCommReq` | `WaitAllSends` / `WaitAllReceives` |

### Pourquoi le filtre marker est appliqué à la réception

Filtrer à l'émission — n'envoyer que les faces du marker clippé — casse la
symétrie entre `sep` et `ext` sur laquelle repose l'appariement par position, et
impose deux mécanismes supplémentaires : transmettre l'identité de chaque face,
et envoyer des messages vides aux pairs que le filtre laisse sans rien à dire
(faute de quoi leur réception bloque).

Le gain est nul — quelques kilo-octets sur un pipeline de 1,5 min — pour deux
invariants fragiles de plus. **Toutes** les faces séparées sont donc envoyées, et
le filtre s'applique à l'arrivée :

```cpp
if(marker == markerBoundaryClipped)
  bdryFaces.emplace_back(*ebf._faceFSDM, marker, fsdmFaceGlobalNumber);
```

---

## Tests

`test/FSClipping/Parallel/FSClippingTestBoundaryExtractionPar.cpp`, enregistré à
**4 procs** — le plus petit compte auquel RCB sépare réellement une paire sur
`cube_hexa_coarse_par.grid`.

| test | vérifie |
|---|---|
| `FaceCountIsPartitionIndependent` | le total vaut 16, la référence séquentielle |
| `NoFaceIsExtractedTwice`          | les clés géométriques sont uniques et totalisent 16 |

Le test a été validé **dans les deux sens** : il échoue (15 ≠ 16) quand l'appel à
`ExchangeSeparatedBdryFaces()` est neutralisé.

> Une première rédaction de ce test passait avec *et* sans le correctif : elle
> comparait le total à une quantité qui variait de la même façon. D'où la
> référence en dur (16) plutôt qu'une valeur recalculée.

### Résultats

| vérification | résultat |
|---|---|
| suite complète | **52/52** |
| balayage 4 markers × 6 comptes de procs (1, 2, 4, 5, 8, 16) | invariant : 16 |
| `mesh_cylinder`, mêmes comptes | invariant : 551 |
| motif cyclique complet (chaque proc → tous), 2 à 16 procs | pas d'interblocage |

---

## Limites connues

- **Le deadlock n'est pas couvert par la suite de tests.** Les maillages du dépôt
  ne produisent aucun échange bidirectionnel (mesuré : 0 paire à 8, 16, 24 et
  32 procs) : ils sont trop petits pour qu'un proc soit à la fois émetteur et
  récepteur vis-à-vis du même voisin. La correction a été validée sur un test
  synthétique du motif d'échange, retiré ensuite. Seul un cas de plusieurs
  millions d'éléments l'exerce réellement.

- **Coût de l'échange non mesuré.** Il porte sur quelques centaines de faces dans
  un pipeline de l'ordre de la minute, donc supposé négligeable — mais aucun
  chiffre ne l'atteste.

- **Invariant à préserver** : les deux moitiés doivent rester triées sur la même
  clé (`VolumeSidePairKey`). Utiliser l'`operator<` de FSDM à la place casse
  l'appariement, comme décrit plus haut.

---

## Fichiers

| fichier | rôle |
|---|---|
| [`src/FSFaceSeparator.cpp`](../src/FSFaceSeparator.cpp) | extraction et échange |
| [`test/FSClipping/Parallel/FSClippingTestBoundaryExtractionPar.cpp`](../test/FSClipping/Parallel/FSClippingTestBoundaryExtractionPar.cpp) | test de non-régression |
| [`test/CMakeLists.txt`](../test/CMakeLists.txt) | enregistrement à 4 procs |

Référence : `negev/infrastructure/src/Distributed/FaceBasedMeshAdapter.cpp`,
méthodes `SeparateFaces()` et `ExchangeSeparatedBdryFaces()`.
