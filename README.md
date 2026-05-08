# 16S rRNA Analysis Pipeline

Snakemake-пайплайн обработки данных 16S рРНК-ампликонного секвенирования: от сырых парных Illumina-ридов до таблицы ASV с таксономией, готовой для анализа разнообразия.

---

## Требования

| Инструмент | Минимальная версия | Примечание |
|---|---|---|
| Snakemake | >= 7.0 | менеджер пайплайна |
| cutadapt | **>= 3.0** | `-j` и `file:` пулы праймеров требуют >= 3.0 |
| Python | >= 3.8 | зависимость Snakemake и cutadapt |
| R + DADA2 | >= 1.26 | через conda-окружение `envs/dada2.yaml` |

```bash
# cutadapt
conda install -c bioconda "cutadapt>=3.0"
cutadapt --version   # должна быть 3.0 или выше

# Snakemake
conda install -c bioconda -c conda-forge snakemake
```

---

## Структура проекта

```
16s/
│
├── Snakefile
├── run_cutadapt.sh
├── README.md
│
├── primers/
│   ├── fwd_primers.fasta
│   └── rev_primers.fasta
│
├── references/                        # gitignored
│   └── silva*.fa.gz
│
├── data/
│   ├── raw/                           # сырые FASTQ (read-only)
│   ├── trimmed/                       # после cutadapt
│   └── qc/                            # после filterAndTrim
│
├── results/
│   ├── derep/                         # дереплицированные объекты (.rds)
│   ├── dada/                          # DADA2 денойзинг (.rds)
│   ├── merged/                        # слитые риды (.rds)
│   ├── filter_stats/                  # статистика фильтрации (.rds)
│   ├── plots/                         # графики моделей ошибок
│   ├── seqtab.rds
│   ├── seqtab_nochim.rds
│   ├── taxa/
│   │   └── taxa.rds
│   └── track.csv                      # учёт ридов по этапам
│
├── logs/
│   ├── cutadapt/
│   ├── dereplication/
│   └── dada_inference/
│
├── scripts/
│   ├── revcomp_primers.py             # обратные комплементы праймеров
│   ├── trim_primers.sh                # вызов cutadapt
│   ├── quality_per_sample.R
│   ├── quality_aggregated.R
│   ├── filter_reads.R
│   ├── error_correction_model.R
│   ├── dereplication.R
│   ├── dada_inference.R
│   ├── merge_reads.R
│   ├── make_seqtable.R
│   ├── remove_chimera.R
│   ├── create_summary.R
│   ├── create_taxa.R
│   └── add_species.R
│
└── envs/
    └── dada2.yaml
```

---

## Установка

```bash
# 1. Клонируй репозиторий
git clone <repo_url>
cd 16s

# 2. Создай структуру папок
mkdir -p data/{raw,trimmed,qc} results/{derep,dada,merged,filter_stats,plots,taxa} \
         logs/{cutadapt,dereplication,dada_inference} references primers

# 3. Положи сырые FASTQ в data/raw/
#    Имена должны соответствовать паттерну {sample}_R1.fastq.gz / {sample}_R2.fastq.gz

# 4. Скачай референсные базы SILVA
wget -P references/ https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
wget -P references/ https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz

# 5. Создай FASTA с праймерами в primers/
```

### Пример `primers/fwd_primers.fasta` (V3–V4)

```
>FWD_341F
CCTACGGGNGGCWGCAG
```

### Пример `primers/rev_primers.fasta` (V3–V4)

```
>REV_805R
GACTACHVGGGTATCTAATCC
```

> Для мультиплексных протоколов добавь все варианты праймеров в эти же файлы — каждый со своим заголовком `>`.

---

## Запуск

### Пробный запуск (без выполнения)

```bash
snakemake --cores 8 --use-conda -n
```

### Полный пайплайн

```bash
snakemake --cores 8 --use-conda
```

### Отдельные шаги

```bash
# Только обрезка праймеров
snakemake --cores 8 --use-conda data/trimmed/sampleA_R1.fastq.gz

# Только графики качества
snakemake --cores 4 --use-conda data/qc/per_sample_forward.pdf

# До удаления химер
snakemake --cores 8 --use-conda results/seqtab_nochim.rds
```

### Параметры обрезки праймеров (`Snakefile → trim_primers`)

| Параметр | По умолчанию | Назначение |
|---|---|---|
| `threads` | `4` | потоков на один cutadapt-вызов |
| `min_len` | `50` | минимальная длина рида после обрезки |
| `error_rate` | `0.1` | допустимая доля ошибок в матче праймера |

Изменяются в блоке `params:` правила `trim_primers` в `Snakefile`.

---

## Этапы пайплайна

| # | Правило Snakemake | Скрипт | Что делает |
|---|---|---|---|
| 1 | `revcomp_primers` | `revcomp_primers.py` | Генерирует обратные комплементы праймеров для отлова read-through |
| 2 | `trim_primers` | `trim_primers.sh` | Обрезает праймеры с обоих концов, выбрасывает необрезанные риды |
| 3 | `quality_per_sample` | `quality_per_sample.R` | Графики качества по каждому образцу |
| 4 | `quality_aggregated` | `quality_aggregated.R` | Агрегированные графики качества |
| 5 | `filter_reads` | `filter_reads.R` | `filterAndTrim`: обрезка по длине и фильтрация по maxEE |
| 6 | `error_correction` | `error_correction_model.R` | `learnErrors`: модель ошибок секвенирования для R1 и R2 |
| 7 | `dereplication` | `dereplication.R` | `derepFastq`: схлопывание одинаковых ридов |
| 8 | `dada_inference` | `dada_inference.R` | `dada`: восстановление ASV |
| 9 | `merge_reads` | `merge_reads.R` | `mergePairs`: слияние R1 и R2 по оверлапу |
| 10 | `make_seqtable` | `make_seqtable.R` | `makeSequenceTable`: матрица образцы × ASV |
| 11 | `remove_chimera` | `remove_chimera.R` | `removeBimeraDenovo`: удаление химер |
| 12 | `assign_taxonomy` | `create_taxa.R` | Наивный байесовский классификатор против SILVA |
| 13 | `add_species` | `add_species.R` | Видовой уровень по точному совпадению |
| 14 | `create_summary` | `create_summary.R` | Таблица учёта ридов по всем этапам → `results/track.csv` |

---

## Ключевые параметры

| Параметр | Где | Что подбирать |
|---|---|---|
| `truncLen = c(240, 200)` | `filter_reads.R` | по графикам качества (где Q падает ниже 30) |
| `maxEE = c(2, 2)` | `filter_reads.R` | строже = меньше ридов, чище данные |
| `minOverlap = 12` | `merge_reads.R` | минимум перекрытия для слияния |
| `threads` | `Snakefile → trim_primers` | потоков cutadapt (по умолчанию 4) |
| `error_rate` | `Snakefile → trim_primers` | допуск ошибок праймера (по умолчанию 0.1) |
