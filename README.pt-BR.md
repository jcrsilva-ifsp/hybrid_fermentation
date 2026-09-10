# hybrid-fermentation-metagenomics

Pipeline de metagenômica shotgun *genome-resolved* para sistemas de fermentação espontânea/mista — recupera MAGs **bacterianos e fúngicos**, anota funcionalmente, e identifica espécies fúngicas por ANI genômico em vez de homologia de proteína.

Desenvolvido por Jean C. R. da Silva no CCBL (FCFRP-USP), sob supervisão do Prof. Ricardo Roberto da Silva.

[Read in English](README.md)

## Motivação

Pipelines padrão de binning metagenômico (ex. MetaWRAP com CONCOCT/MaxBin2/MetaBAT2) foram desenhados quase exclusivamente para procariontes. Genomas fúngicos — maiores, com íntrons e uso de códons distinto — tendem a ser fragmentados ou descartados no binning. Em bebidas onde a fração fúngica é central (leveduras em Pulque, Kombucha, Lambic), isso produz um ponto cego sistemático: o pipeline nunca recupera um único MAG fúngico, mesmo com DNA presente em abundância.

Este repositório implementa e valida uma correção: bifurcação de domínios (Tiara) antes do binning, para que bactérias e fungos sejam processados por ferramentas apropriadas em ramos separados. Validado contra três datasets públicos independentes (Pulque, Kombucha, Lambic — ver abaixo).

Um segundo problema, independente do primeiro, foi encontrado e corrigido durante a validação: identificação taxonômica de MAGs fúngicos por homologia de proteína contra o Swiss-Prot identifica errado leveduras de fermentação sub-representadas nesse banco (ex. *Brettanomyces bruxellensis* praticamente não tem entradas de genoma nuclear no Swiss-Prot). Comparação de ANI genômico (fastANI) contra um pequeno conjunto de genomas de referência resolve isso, sem depender da cobertura de anotação do banco de referência.

## Visão geral do pipeline

| Fase | Ferramenta(s) | Saída |
|---|---|---|
| 1 | Trimmomatic (+ bmtagger opcional, só paired-end) | Reads limpas |
| 1c | Kraken2 + Bracken | Perfil taxonômico direto |
| 2 | MEGAHIT | Contigs montados |
| 3 | Tiara | Contigs bifurcados (procariótico/eucariótico) |
| 4A | MetaWRAP (CONCOCT+MaxBin2+MetaBAT2) + CheckM | MAGs bacterianos |
| 4B | MetaBAT2 + BUSCO (fungi_odb12.2) | MAGs fúngicos |
| 5 | Prodigal + PROKKA (bactérias); MetaEuk (fungos); Kofamscan/KEGG (ambos) | Anotação funcional |
| 5b | fastANI contra conjunto de referências | Identificação de espécie fúngica |
| 6 | antiSMASH + dbCAN | Prospecção de BGCs/CAZymes (não configurada nesta versão) |

Duas variantes: `pipeline_hibrido_corrigido.sh` (paired-end, Illumina) e `pipeline_hibrido_singleend.sh` (single-end, ex. Ion Torrent).

## Instalação

Seis ambientes conda/micromamba isolados são necessários (ver `environments/*.yml`).

```bash
git clone <REPO_URL>
cd hybrid-fermentation-metagenomics

for ENV in tiara metaeuk busco kofamscan bracken fastani; do
    micromamba env create -p ./env_${ENV} -f environments/env_${ENV}.yml
done
```

O MetaWRAP **não** é instalável via conda de forma confiável (o nome do pacote `metawrap` colide com um pacote Python não relacionado). Instale a partir do repositório oficial:

```bash
git clone https://github.com/bxlab/metaWRAP.git
# siga as instruções de instalação do próprio repositório (config-metawrap)
```

### Dados de referência

| Dado | Uso | Origem |
|---|---|---|
| Swiss-Prot (indexado MMseqs2/MetaEuk) | Predição gênica fúngica por homologia | UniProt |
| fungi_odb12.2 | Completude de MAGs fúngicos (BUSCO) | Downloads oficiais BUSCO |
| Perfis KOfam + ko_list | Mapeamento funcional KEGG | ftp.genome.jp/pub/db/kofam/ |
| Kraken2 PlusPF | Perfil taxonômico direto | Coleção padrão Kraken2 |
| Genomas de referência fúngicos para ANI (Fase 5b) | Identificação de espécie fúngica | NCBI (baixados automaticamente — ver abaixo) |

A Fase 5b baixa os genomas de referência automaticamente via `datasets` (NCBI Datasets CLI, incluído no `env_fastani`) na primeira execução, usando os *accessions* configurados no topo do script:

```bash
ANI_REFS["Saccharomyces_cerevisiae"]="GCF_000146045.2"
ANI_REFS["Brettanomyces_bruxellensis"]="GCF_011074885.1"
ANI_REFS["Zygosaccharomyces_rouxii"]="GCA_000026365.1"
```

Adicione ou remova entradas conforme as leveduras relevantes ao seu sistema (ex. *Pichia*, *Hanseniaspora*, *Dekkera*, cepas adaptadas a lager).

### Reproduzindo o ambiente exato usado neste trabalho

```bash
micromamba activate -p ./env_<nome>
conda env export | grep -v "^prefix: " > environments/env_<nome>.yml
```

## Uso

Ajuste `BASE_DIR` e cada placeholder `/path/to/...` no topo do script para o seu ambiente. Depois, submeta um job por amostra:

```bash
sbatch --export=ALL,SAMPLE_ID=<accession> --job-name=metag_<accession> pipeline_hibrido_corrigido.sh
```

**Fortemente recomendado**: rode uma amostra piloto primeiro, acompanhando o log (`$BASE_DIR/logs/<accession>.log`), antes de submeter em lote. Diretivas de recurso do SLURM (partição, CPUs, memória) são específicas do cluster e precisam ser reconferidas em qualquer ambiente novo.

## Validação

Testado contra três datasets públicos independentes de fermentação espontânea/mista:

- **Pulque** (Chacón-Vargas et al., 2020, *Scientific Reports*, BioProject PRJNA603591): taxonomia de MAGs bacterianos consistente com o artigo de referência; MAGs fúngicos recuperados e identificados corretamente como *Saccharomyces cerevisiae* (ANI ~87-88% contra a referência S288C — abaixo do limiar convencional de ≥95% para espécie, provavelmente refletindo divergência real de linhagem, não incerteza de identificação); *Zymomonas mobilis*, ausente dos MAGs, confirmado presente e com dinâmica consistente com o artigo via Kraken2+Bracken.
- **Kombucha** (Suhre et al., 2025, *Food Bioscience*): MAGs fúngicos, originalmente identificados errado como *S. cerevisiae* por homologia de proteína, corrigidos para *Brettanomyces bruxellensis* por ANI (muitos hits ≥94-98%) — consistente com o artigo de referência, que relata *Brettanomyces* dominando a comunidade fermentada (>88% em algumas amostras).
- **Lambic** (De Roos et al., 2020, *Frontiers in Microbiology*, ENA BioProject PRJEB28363; Ion Torrent single-end): sucessão bacteriana (tipo-Acetobacter → tipo-Lactobacillales ao longo do tempo de fermentação) reproduz a dinâmica de acidificação seguida de maturação descrita no artigo; identificação fúngica por ANI: 3/6 amostras confirmadas *S. cerevisiae* (ANI ~99%), 1/6 corrigida para *Brettanomyces bruxellensis* (ANI 98,45%), 2/6 não resolvidas contra o conjunto de referências atual (candidata: *Pichia membranifaciens*, relatada no artigo mas ainda não incluída como referência).

## Limitações conhecidas

- `metawrap bin_refinement` falha com frequência (bug da própria ferramenta); um fallback via CheckM lida com isso automaticamente.
- Nem toda amostra produz MAG fúngico — ausência é um resultado válido, não falha.
- O conjunto de referências da Fase 5b é necessariamente incompleto; um MAG sem correspondência confiável de ANI pode só precisar de um genoma de referência ainda não incluído.
- Fase 6 (antiSMASH + dbCAN) não configurada nesta versão.
- O Kraken2 PlusPF não inclui genomas de plantas; em amostras com matriz vegetal (ex. Kombucha com fruta), grande parte das reads pode ficar "não classificada" por esse motivo, não por falha do método.

## Estrutura do repositório

```
.
├── README.md
├── README.pt-BR.md
├── pipeline_hibrido_corrigido.sh       (paired-end)
├── pipeline_hibrido_singleend.sh       (single-end)
├── environments/
│   ├── env_tiara.yml
│   ├── env_metaeuk.yml
│   ├── env_busco.yml
│   ├── env_kofamscan.yml
│   ├── env_bracken.yml
│   └── env_fastani.yml
└── LICENSE
```

## Autoria

Jean Carlos Rodrigues da Silva, pós-doutorando — CCBL, FCFRP-USP.
Supervisão: Prof. Ricardo Roberto da Silva.

## Licença

A definir (ver `LICENSE`).
