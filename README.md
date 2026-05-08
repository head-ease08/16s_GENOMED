#  16S rRNA Analysis Pipeline

Полный пайплайн обработки данных 16S рРНК-ампликонного секвенирования: от сырых парных Illumina-ридов до `phyloseq`-объекта с филогенетическим деревом, готового для анализа разнообразия и таксономического состава.

---

## 📁 Структура проекта

```
16short/
│
├── 📄 README.md
├── 📑 dada2_paired_cutadapt_tree.Rmd     # основной R-пайплайн
├── 🐚 run_cutadapt.sh                    # обрезка праймеров
│
├── 🧬 primers/
│   ├── fwd_primers.fasta
│   └── rev_primers.fasta
│
├── 📋 metadata/
│   ├── sample-metadata.tsv
│   └── manifest.tsv
│
├── 📚 references/                        # gitignored
│   ├── silva_nr_v132_train_set.fa.gz
│   └── silva_species_assignment_v132.fa.gz
│
├── 💾 data/
│   ├── raw/                              # сырые FASTQ (read-only)
│   ├── trimmed/                          # после cutadapt
│   │   ├── logs/
│   │   └── cutadapt_summary.tsv
│   └── qc/                               # после filterAndTrim
│
├── 📈 results/
│   ├── track.csv                         # учёт ридов по этапам
│   ├── seqtab_nochim.rds
│   ├── taxa.rds
│   ├── asv_tree.nwk
│   ├── ps.rds                            # финальный phyloseq
│   └── plots/
│
└── 🌐 MicrobiomeAnalyst_data/
    ├── asv_table.csv
    ├── tax_table.csv
    ├── rep_seqs.csv
    ├── metadata.csv
    └── asv_tree.nwk
```

---

## Requirements

| Tool | Minimum version | Notes |
|------|----------------|-------|
| cutadapt | **>= 3.0** | `-j` multithreading and `file:` adapter pools require >= 3.0 |
| Python | >= 3.6 | cutadapt runtime dependency |

Install cutadapt:

```bash
pip install "cutadapt>=3.0"
# or via conda:
conda install -c bioconda "cutadapt>=3.0"
```

Verify before running:

```bash
cutadapt --version   # must print 3.0 or higher
```

---

## 📥 Установка

```bash
# 1. Клонируй репозиторий
git clone <repo_url> dead_mans_teeth
cd dead_mans_teeth

# 2. Создай структуру папок
mkdir -p data/{raw,trimmed,qc} results/plots references primers metadata

# 3. Положи сырые FASTQ в data/raw/
#    Имена должны соответствовать паттерну *_R1*.fastq(.gz) / *_R2*.fastq(.gz)

# 4. Скачай референсные базы SILVA
wget -P references/ https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
wget -P references/ https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz

# 5. Положи метаданные в metadata/
#    sample-metadata.tsv — основная таблица образцов
#    manifest.tsv        — соответствие имён файлов и образцов

# 6. Создай FASTA с праймерами в primers/
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

> 💡 Для мультиплексных протоколов добавь все варианты праймеров в эти же файлы — каждый со своим заголовком `>`.

---

## 🚀 Запуск

### Шаг 1: обрезка праймеров

```bash
chmod +x run_cutadapt.sh

./run_cutadapt.sh \
    data/raw/ \
    data/trimmed/ \
    primers/fwd_primers.fasta \
    primers/rev_primers.fasta
```

```bash
cat data/trimmed/cutadapt_summary.tsv
```

---

## 📚 Описание этапов

### 1️⃣ Cutadapt — удаление праймеров

Скрипт `run_cutadapt.sh` запускает cutadapt в режиме **linked adapters**:

- 🔍 ищет forward-праймер в начале R1 и reverse-праймер в начале R2
- 🔄 дополнительно ловит read-through (обратные комплементы на 3'-концах)
- 🗑️ выкидывает риды без праймера через `--discard-untrimmed`
- 📝 логирует каждый образец в `data/trimmed/logs/<sample>.log`

### 2️⃣ Quality filtering

`filterAndTrim` обрезает риды по длине (`truncLen`) и фильтрует по `maxEE` (ожидаемые ошибки). Параметры выбираются по графикам качества.

### 3️⃣ Error model

`learnErrors` обучает параметрическую модель ошибок секвенирования отдельно для R1 и R2.

### 4️⃣ DADA2 inference

Алгоритм восстанавливает истинные биологические последовательности (ASV) с разрешением до одного нуклеотида.

### 5️⃣ Merge pairs

`mergePairs` склеивает forward и reverse ASV по оверлапу (минимум 12 п.н., 0 mismatches по умолчанию).

### 6️⃣ Chimera removal

`removeBimeraDenovo` с методом `consensus`. Доля химер обычно 5–20%.

### 7️⃣ Taxonomy

Наивный байесовский классификатор против SILVA. Видовой уровень ненадёжен для 16S — используй как ориентир.

---

## ⚙️ Параметры и тюнинг

### Переменные окружения для `run_cutadapt.sh`

|Переменная|Дефолт|Назначение|
|---|---|---|
|`THREADS`|`4`|потоков на один cutadapt-вызов|
|`MIN_LEN`|`50`|минимальная длина рида после обрезки|
|`ERROR_RATE`|`0.1`|допустимая доля ошибок в матче праймера|

Пример:

```bash
THREADS=8 MIN_LEN=100 ./run_cutadapt.sh data/raw/ data/trimmed/ \
    primers/fwd_primers.fasta primers/rev_primers.fasta
```

### Ключевые параметры в Rmd

|Параметр|Где|Что подбирать|
|---|---|---|
|`truncLen = c(240, 200)`|`filterAndTrim`|по графикам качества (Q ниже 30)|
|`maxEE = c(2, 2)`|`filterAndTrim`|строже = меньше ридов, чище данные|
|`minOverlap = 12`|`mergePairs`|минимум перекрытия для склейки|
|`seqtab.lenfilt` диапазон|после `mergePairs`|под ожидаемый размер амплкона|
