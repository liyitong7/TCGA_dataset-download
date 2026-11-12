# TCGA_dataset-download
TCGA数据集下载器
项目来源：https://github.com/NCI-GDC/gdc-client
不太喜欢带参执行的用法，复刻后省略了固定参数，仅限下载TCGA数据集用。

1.打开GDC Data Potral： https://portal.gdc.cancer.gov/analysis_page?app=Downloads
2.点击页面中间的"Projects"，或者直接在搜索框搜 TCGA
3.根据你的项目需求(做泛癌分析)，你可能需要勾选多个项目(比如 TCGA-BRCA,TCGA-LUAD 等等)，或者干脆选整个 TCGA 项目。
4.选好后，点击"View Files in Repository"(查看仓库中的文件)。
5.在左侧的过滤条件(Filters)里:

    DataFormat选tsv或txt(通常表达量数据是 tsv)。
    Data Category选Transcriptome Profiling 。
    DataType选Gene Expression Quantification 。
    ExperimentalStrategy选RNA-Seq 。
    WorkflowType选 STAR-counts(这是目前最主流的)

6.过滤完以后，页面会显示几千上万个文件。（Total of 11505 Files， 10517 Cases， 48.74 GB）
点击页面右上角全选(或者按需勾选)，然后点击"Add toCart"(加入购物车)。
或者点击Mainfest，下载好文件清单，然后进行下载。
