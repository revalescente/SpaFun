# SpaFun — progetto di sviluppo

*Bozza dell'8 ottobre 2026, scritta a partire dal ramo `origin/dev` (commit `7f58c46`).*

**Obiettivo:** un pacchetto Bioconductor generale per l'analisi di pattern di punti spaziali multitipo (cellule) con metodi di analisi di dati funzionali (FDA). Deve funzionare su qualsiasi tecnologia di spatial omics a risoluzione di singola cellula (H&E segmentato, CODEX/IMC/MIBI, Xenium/CosMx, Visium HD), con TCGA come caso d'uso guida e la scala dei vetrini interi (milioni di cellule) come punto di forza.

---

## 1. Dove siamo

Nel ramo `dev`:

| Livello | Funzioni | Stato |
|---|---|---|
| Lettura | `read_cells_h5ad()`, `cells2ppp()`, `ann2ppp()` | Rivisto, testato, verificato su TCGA |
| Downsampling | `downsample()` + `_uniform` / `_KDE` / `_window` | Rivisto, testato |
| Summary functions | — (solo codice nelle vignette: `Lcross` scritto a mano) | Da fare |
| FDA | `fsmooth()`, `plot_smoothed_Lfun()` | Codice originale, non rivisto, senza test |
| Analisi | fPCA e clustering nelle vignette | Codice delle vignette, nessuna API |
| Dati | `tcga_crops`, `inst/extdata/tcga_tiny.*` | OK (sotto i 5 MB) |
| Infrastruttura | fusen, pkgdown, GitHub Actions, testthat | Presente, con lacune (§2.4) |

Nota: `main` è indietro di due commit rispetto a `dev`.

---

## 2. Altre sistemazioni (prima di aggiungere funzionalità)

In ordine di priorità.

### 2.1 Errori e problemi di correttezza

1. **`fsmooth()` funziona solo con data.table.** `dt[, -1, with = FALSE]` fallisce su un data.frame, e `data.table` sta solo in Suggests. Basta usare `as.matrix(dt[, -1])` (o `dt[, -1, drop = FALSE]`) e accettare anche una matrice con `r` a parte.
2. **Unità di misura.** Le coordinate sono in pixel a 40x, ma `fsmooth()` scrive "micron" nei `fdnames`. A circa 0,25 µm/pixel, `r = 0:150` vale circa 0–37 µm, non 0–150. Serve un'unica scelta esplicita:
   - un argomento `mpp` (microns per pixel) oppure `unit_scale` in `cells2ppp()`, e `unitname(pp) <- c("micron", "microns")`;
   - poi convertire tutto in µm, così `r` ha un significato biologico ed è confrontabile tra vetrini e tecnologie diversi.
3. **Fallback silenzioso in `fsmooth()`.** Se `smooth.basis` fallisce, il modello passa a `Lfdobj = 2` e a un `lambda` più grande, ma il messaggio compare solo con `verbose`. Deve diventare un `warning()`, e i parametri usati davvero vanno salvati nel risultato.
4. **`lambda` fisso a 1e-6.** Va scelto per GCV (`lambda = "gcv"` come default, con una griglia su `log10(lambda)`), altrimenti il grado di smoothing dipende dalla scala delle curve.
5. **`plot_smoothed_Lfun()`:**
   - se `clusters_vec` non ha nomi, tutti i cluster diventano `NA` senza avvisi;
   - con `fsmooth(centered = TRUE)` e `center_plot = TRUE` la centratura viene fatta due volte (il risultato di `fsmooth` deve ricordare se è già centrato);
   - il tag roxygen è rotto (`#' #' @importFrom utils globalVariables`);
   - "L-cross" è scritto fisso nel titolo.
6. **`fdnames` in italiano** ("Distanza", "Funzione") dentro il pacchetto: vanno in inglese o presi dai metadati.
7. **Tipo cellulare mancante.** Se un tipo non c'è in un vetrino, `Lcross` va in errore o produce NA. Oggi lo gestisce solo il codice delle vignette: deve gestirlo la funzione (§3.3).

### 2.2 API e naming (adesso cambiare non costa nulla, dopo il rilascio sì)

- Usare un solo stile, snake_case minuscolo: `downsample_kde()`, `plot_fd()`/`plot_curves()`, `smooth_curves()` (al posto di `fsmooth`), `h5ad_to_ppp()` o semplicemente `read_cells()`. Il nome `ann2ppp()` va bene come alias deprecato.
- Argomenti con nomi descrittivi: `M` → `norder`, `genLfun.fd` → `Lfdobj` o `penalty_order`.
- Risultati con una classe S3 (`spafun_fd`, con metodi `print`, `plot`, `as.data.frame`) invece di liste semplici.
- Ogni funzione esportata deve avere un `@examples` eseguibile sui dati di esempio: BiocCheck lo richiede.

### 2.3 Metadati del pacchetto

- `Description:` ora dice "Function to analyse TCGA dataset". Va riscritta in 3–5 frasi sul metodo.
- `Title:` va in Title Case e senza abbreviazioni ("Functional Data Analysis of Multitype Spatial Point Patterns").
- `biocViews:` è vuoto. Proposta: `Software, Spatial, SingleCell, StatisticalMethod, Clustering, Regression, Visualization, ImmunoOncology`.
- `Depends: R (>= 3.5)` → la versione R della Bioc devel al momento dell'invio.
- Email: in DESCRIPTION c'è `@phd.unipd.it`, l'account usa `@unipd.it`. Scegline una stabile nel tempo, perché il maintainer Bioc la usa per anni.
- `rhdf5` e `data.table` in Suggests: vanno tenuti solo se usati in vignette o test.
- Provenienza e licenza dei dati di esempio (segmentazioni HoVer-Net dei TCGA, da quale fonte?) vanno documentate in `man/tcga_crops.Rd` e `data-raw/`.

### 2.4 Infrastruttura

- `R-CMD-check.yaml` parte solo su `main`/`master`: i push su `dev` non vengono controllati. Aggiungi `dev`, oppure lavora con PR verso `main`.
- Aggiungi un job con **BiocCheck** (`bioc-actions` o `r-lib/actions` + `BiocCheck::BiocCheck()`).
- Mancano test per `fsmooth()` e `plot_smoothed_Lfun()`. Obiettivo: copertura sopra l'80% (codecov è già configurato, `codecov: false` in fusen).
- Il README è ancora il template (`summary(cars)`): serve un esempio reale di 10 righe e il badge della CI.
- Unisci `dev` → `main` e assegna un tag `v0.1.0`.
- **fusen**: va bene anche per Bioconductor, a patto che `R/`, `tests/` e `vignettes/` siano sempre gonfiati e allineati prima di ogni push. Un check in CI che rigonfia e verifica `git diff --exit-code` lo garantisce.

### 2.5 Pipeline nelle vignette

- I percorsi HPC (`/projects/shared/TCGA_data/...`) vanno tolti dalle vignette, anche con `eval: false`. Il workflow HPC va spostato nel repo `ppp_TCGA` (meglio ancora con `targets` + `crew.cluster`/`future.batchtools`); il pacchetto deve offrire solo funzioni.
- Il codice di fPCA, clustering e grafici dei punteggi (`par(mfrow)`, `plot` base) va trasformato in funzioni (§3.4).

---

## 3. Architettura proposta

```
 Input                Preprocessing          Summary functions        Rappresentazione        Analisi
 ─────                ─────────────          ─────────────────        funzionale              ──────
 SpatialExperiment ─┐  finestra/maschera      K, L, g, J, mark-        smoothing (GCV)         descrittiva: fPCA,
 SingleCellExp.    ─┼→ unità (µm)          →  connection, cross     → registrazione        →  MFPCA, clustering
 AnnData (h5ad)    ─┤  QC (tipi rari, NA)     tutte le coppie i,j      fd multivariato         inferenziale: FoSR,
 data.frame / ppp  ─┘  downsampling           inhom., inviluppi        classe dedicata         SoFR/Cox, test, LMM
                                              ↓                                                 ↓
                                     SummarizedExperiment (r × campioni)              grafici (ggplot2)
```

### 3.1 Modello dati (decisivo per Bioconductor)

I revisori Bioc chiedono di **riusare le classi esistenti** invece di inventarne di nuove.

- **Input:** `SpatialExperiment` (`spatialCoords()`, tipo cellulare in `colData`, `sample_id` per i vetrini). È la porta d'ingresso naturale per lavorare su più campioni. `anndataR` converte già h5ad → `SingleCellExperiment`; l'attuale `read_cells_h5ad()` resta come lettore veloce.
- **Summary functions stimate:** un `SummarizedExperiment` con
  - righe = combinazioni (funzione, tipo_i, tipo_j, r), `rowData` con queste colonne;
  - colonne = campioni/vetrini, `colData` = dati clinici (sopravvivenza, stadio, ecc.);
  - `assays`: `estimate`, eventualmente `theo` (CSR), `lo`/`hi` (inviluppi), `n_i`, `n_j`.

  Così i dati clinici viaggiano insieme alle curve, il subsetting `se[, se$stage == "III"]` funziona subito e l'oggetto si integra con il resto dell'ecosistema.
- **Oggetto funzionale:** una classe S3 leggera sopra `fda::fd` (o `funData`, se adotti MFPCA), con un riferimento ai metadati.

### 3.2 Input, preprocessing, scala

- `spe_to_ppp(spe, sample_id, type_col, window = c("mask", "convex", "concave", "bbox"))`. Per le tecnologie senza maschera serve una finestra stimata: convex hull, alpha-shape (`concaveman`) o la maschera dal tessuto.
- Unità: §2.1.2.
- QC: soglia minima di cellule per tipo (`min_n`), con un rapporto di cosa è stato escluso, nello stile di `n_removed`.
- Parallelismo: oggi usi `future.apply`. Bioc preferisce **BiocParallel** (`BPPARAM`) e alcuni revisori lo chiedono esplicitamente. Valuta di passare a `BiocParallel::bplapply`: con `BatchtoolsParam` funziona anche su SLURM.

### 3.3 Motore delle summary functions (nuovo, priorità alta)

```r
sumfun(x, fun = c("L", "K", "g", "J", "markconnect"),
       types = NULL,            # NULL = tutte le coppie
       pairs = c("cross", "all", "self"),
       r = NULL, correction = "border",
       inhom = FALSE,           # Lcross.inhom con intensità stimata
       downsample = list(method = "uniform", n_target = 4e4),
       envelope = c("none", "csr", "random_labelling"),
       BPPARAM = SerialParam())
# -> SummarizedExperiment
```

- Griglia `r` comune e condivisa, scelta in µm (default: fino a ¼ del lato più corto della finestra, come `rmax.rule`).
- Coppie mancanti → NA con un flag, invece di un errore.
- **Random labelling** come riferimento per le funzioni cross. In biologia è più sensato di CSR, perché la disposizione spaziale complessiva delle cellule non è di interesse. Restituire `L_ij − L_ij^{RL}` elimina l'effetto dell'architettura del tessuto. Probabilmente è il contributo metodologico più utile da avere già pronto.
- **Inomogeneità:** con vetrini interi l'intensità non è costante, quindi va offerta `Lcross.inhom` / `Kinhom` con l'intensità stimata (anche qui il downsampling uniforme è compatibile).
- Prestazioni: per milioni di cellule, `downsample_window` (tile) più la media delle stime per tile è un'alternativa al thinning che dà anche una stima della varianza.

### 3.4 Livello FDA descrittivo

- `smooth_curves(se, ...)`: B-spline con GCV; opzione di vincolo `L(0) = 0`.
- `fpca(fdobj, nharm = NULL, pve = 0.9)`, con scelta del numero di componenti per quota di varianza spiegata.
- **Multivariata (il cuore "alta dimensionalità"):** con K tipi cellulari ci sono K(K−1)/2 curve cross per campione (6 tipi → 15 curve). Opzioni:
  - MFPCA (pacchetto `MFPCA`, Happ & Greven) sulle curve concatenate;
  - fPCA per coppia seguita da PCA/clustering sui punteggi impilati;
  - selezione sparsa delle coppie informative (group lasso sui punteggi).
- Clustering: k-means/Ward sui punteggi (come ora), più metodi funzionali veri (`funHDDC`, `fdakmeans`), con scelta di k (silhouette, gap).
- Registrazione (opzionale): allineamento delle curve se il raggio di interazione varia tra campioni.

### 3.5 Livello FDA inferenziale

- **Test tra gruppi:** ANOVA funzionale con permutazioni (`fdANOVA` o una versione propria), bande simultanee.
- **Regressione scalare su funzione:** outcome clinico ~ curve (`refund::pfr`); **Cox funzionale** per la sopravvivenza TCGA (`refund`/`mgcv` con `family = cox.ph`).
- **Funzione su scalare** (FoSR): come cambia L_ij con stadio/sottotipo.
- **Più ROI o vetrini per paziente:** modelli funzionali misti (`refund::pffr`, oppure `mgcv` con effetti casuali). È frequente in imaging multiplex.
- Le stime di L su punti diradati hanno varianza diversa a seconda di n: valuta pesi per campione nella regressione.

### 3.6 Grafici

`plot_curves()`, `plot_fpca()` (armoniche ± sd, punteggi), `plot_clusters()`, `plot_ppp()` (vetrino con i tipi colorati, a bassa risoluzione grazie al downsampling KDE), `plot_pairs_heatmap()` (riassunto di tutte le coppie i–j). Tutti in ggplot2, senza `par()`.

---

## 4. Posizionamento rispetto ai pacchetti esistenti

Da verificare prima di scrivere introduzione e paper: non ho potuto controllare lo stato attuale online, quindi le informazioni qui sotto vanno riscontrate.

| Pacchetto | Cosa fa | Differenza possibile di SpaFun |
|---|---|---|
| **mxfda** (CRAN, Wrobel e coll.) | Summary functions K/L/g, FPCA, Cox funzionale su imaging multiplex | **Il concorrente più vicino.** Va studiato a fondo. Differenze: scala dei vetrini interi, downsampling con garanzie, multivariata su tutte le coppie, integrazione SpatialExperiment/Bioc |
| **spicyR / lisaClust** (Bioc) | L-function cross + modelli misti; cluster locali | Usano un riassunto scalare della curva, non l'intera funzione |
| **SpaceANOVA** | ANOVA funzionale su K/g multitipo | Solo test di gruppo |
| **scider, SPIAT, imcRtools, spatialTIME** | Statistiche spaziali varie per singola cellula | Nessun approccio FDA sistematico |
| **Voyager / SpatialFeatureExperiment** | Statistica geospaziale (Moran, ecc.) | Prospettiva diversa (lattice/geostatistica, non point process) |

**Messaggio proposto:** *"Analisi FDA completa e scalabile delle summary functions multitipo, dai vetrini interi alla regressione con outcome clinici, nell'ecosistema Bioconductor."*

---

## 5. Roadmap per fasi

Le durate sono indicative e pensate per lavoro part-time accanto al dottorato.

### Fase 0 — Consolidamento (2–3 settimane) → `v0.1.0`
- [ ] Tutto il §2 (errori, naming, metadati, CI su `dev`, BiocCheck in CI, test per FDA e grafici, README)
- [ ] Unità in µm per tutto il pacchetto
- [ ] Merge `dev` → `main`, tag

### Fase 1 — Modello dati Bioc e lettura (3–4 settimane) → `v0.2.0`
- [ ] `SpatialExperiment` in ingresso, `spe_to_ppp()`, finestre stimate
- [ ] `tcga_crops` anche come `SpatialExperiment` di esempio
- [ ] Decisione `future` vs `BiocParallel`

### Fase 2 — Motore delle summary functions (4–6 settimane) → `v0.3.0`
- [ ] `sumfun()` → `SummarizedExperiment`, tutte le coppie, NA gestiti
- [ ] Random labelling e versioni inomogenee
- [ ] Benchmark su HPC: tempo e memoria per 10⁵ / 10⁶ / 5·10⁶ cellule, thinning vs tile
- [ ] Rifare i risultati TCGA con la nuova pipeline (i vecchi erano con KDE)

### Fase 3 — FDA descrittiva e multivariata (4–6 settimane) → `v0.4.0`
- [ ] `smooth_curves()` con GCV, classe `spafun_fd`
- [ ] `fpca()`, MFPCA, clustering con scelta di k
- [ ] Grafici
- [ ] Studio di simulazione: processi multitipo noti (Thomas, Strauss multitipo, pattern segregati) → la pipeline ritrova i gruppi?

### Fase 4 — FDA inferenziale (6–10 settimane) → `v0.5.0`
- [ ] Test funzionali tra gruppi, FoSR, SoFR/Cox
- [ ] Applicazione ai TCGA con la sopravvivenza (materiale per il **paper metodologico**)
- [ ] Validazione su un secondo dataset di tecnologia diversa (es. un dataset IMC/CODEX pubblico già disponibile in `SpatialExperiment` tramite ExperimentHub)

### Fase 5 — Invio a Bioconductor (3–4 settimane) → `0.99.0`
- [ ] Checklist §6 tutta verde
- [ ] Invio per uno dei due cicli annuali di rilascio (aprile/ottobre; controlla la scadenza di invio sul sito Bioc, di solito qualche settimana prima del rilascio)
- [ ] Dopo l'accettazione: **application note** (Bioinformatics / JOSS / F1000Research Bioconductor gateway)

**Ordine consigliato per i paper:** il paper metodologico usa le versioni 0.4–0.5 da GitHub. L'application note arriva dopo l'accettazione Bioc, che fa da revisione del software.

---

## 6. Checklist Bioconductor

- [ ] `R CMD check` senza ERROR/WARNING e con NOTE giustificate; `BiocCheck()` senza ERROR
- [ ] Versione `0.99.0` all'invio; `R (>= <versione devel>)`
- [ ] `biocViews` validi; `Description` ≥ 3 frasi
- [ ] Usa classi Bioc (`SpatialExperiment`/`SummarizedExperiment`) in ingresso e/o in uscita
- [ ] Vignetta con contenuto vero ed eseguibile, con sezione `sessionInfo()`; check totale entro pochi minuti
- [ ] `@examples` eseguibili per ogni funzione esportata (niente `\dontrun` senza motivo)
- [ ] Test unitari (testthat), buona copertura
- [ ] Tarball sotto i 5 MB; dati più grandi → pacchetto dati separato su **ExperimentHub** (es. `SpaFunData` con alcuni vetrini TCGA completi)
- [ ] Stile: righe ≤ 80 caratteri, indentazione 4 spazi, niente `T`/`F`, niente `sapply` (BiocCheck li segnala), funzioni non troppo lunghe
- [ ] Niente `print()`/`cat()` per i messaggi: `message()`/`warning()`
- [ ] `NEWS.md` aggiornato
- [ ] Maintainer iscritto a bioc-devel e al Support Site, con tag del pacchetto
- [ ] Repo GitHub pubblico; invio tramite issue su `Bioconductor/Contributions`
- [ ] Licenza compatibile (GPL-3 va bene); provenienza dei dati di esempio documentata

---

## 7. Decisioni aperte

1. `fda` vs `funData`/`MFPCA` vs `refund` come base FDA: probabilmente `fda` per lo smoothing e `refund` per la regressione, ma la classe interna va decisa una volta sola.
2. `future` vs `BiocParallel`.
3. Mantenere `fusen` a lungo termine o tornare a `R/` classico dopo la v0.1? Per i contributori esterni `R/` è più familiare.
4. Nome: "SpaFun" è libero su CRAN/Bioc? Verificare. Va bene anche per il paper?
5. Pacchetto dati separato (ExperimentHub) sì/no.
