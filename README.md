# 인간 단백질 상호작용(Human PPI) 네트워크 기반 서열 데이터셋 
해당 프로젝트는 인간 단백질 상호작용 네트워크 데이터셋(1)으로 부터 데이터마이닝, 머신러닝 그리고 딥러닝 분야에서 다양하게 활용할 수 있도록 아미노산 서열, PSSM(Position Specific Score Matrix)를 정제하는 작업을 소개하기 위해 제작 되었습니다.

## SNAP - Human Protein-Protein Interaction(PPI) Network 
스탠포드 네트워크 분석 플랫폼([Stanford Network Analysis Platform]([http://snap.stanford.edu/index.html](http://snap.stanford.edu/index.html)))에서는 다양한 네트워크 구조 데이터셋을 제공하고 있습니다.  이 중 BioSNAP 프로젝트는 Network Biology분야에서 활용할 수 있는 네트워크 데이터셋을 제공합니다. 대부분 Network 기반 분석은 노드 사이의 연결(edge)정보를 바탕으로 이루어지므로 각 노드에 대한 자세한 데이터는 제공하지 않습니다. 여기서 사용할 데이터셋은 Human Protein-Protein Interaction Network (PP-pathway)로, 인간 단백질 (gene으로 불리기도 함)의 ID값 쌍으로 네트워크 edge가 구성되어 있습니다. 

    # PP-pathways_ppi.csv
    1394,2778 # pair-1(gene ID, gene ID)
    6331,17999 # pair-2(gene ID, gene ID)
    ...
    4790,79155 # pair-N(gene ID, gene ID)
정수 값의 쌍으로 구성되는 edge list가 바로 SNAP에서 제공하는 네트워크 데이터셋 입니다.  이 정수 값들은 NCBI gene ID이며 entrez ID로 불립니다. 학계에서 발견된 gene에 고유한 식별자를 부여하기 위해 NCBI에서 지정한 값 입니다. 다음 단계에서 우리는 각 gene ID에 대응하는 아미노산 서열 정보를 얻을 것 입니다.  서열 정보는 생물정보학 분야에서 머신러닝 응용 연구에서 많이 사용됩니다. 대표적으로 PPI classification 또는 RNA-Protein Interaction(RPI) Network 에서 RPI classification으로 연구분야가 존재하며, 입력 데이터로 서열 정보를 사용합니다. 해당 프로젝트에서는 PP-pathway를 PPI classification과 같은 머신러닝 연구에 활용할 수 있도록 정제하는 과정 다룰 예정입니다.

## Entrez ID를 UniProt ID로 변경하기
PP-pathway에 존재하는 모든 노드 gene ID의 아미노산 서열을 구하기 위해 [UniProt]([https://www.uniprot.org/](https://www.uniprot.org/)) 데이터베이스를 이용합니다.  [UniProt 다운로드 센터]([https://www.uniprot.org/downloads](https://www.uniprot.org/downloads))에서 [검토가 완료된 Swiss-Prot 단백질 서열 fasta 파일 (Reviewd Swiss-Prot)](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz)을 다운로드 받습니다. Fasta 파일은 ID와 sequence 가 연속적으로 기록된 포맷으로, 앞서 다운로드 받은 서열 데이터를 참고하여 PP-pathway 데이터셋의 노드들의 서열을 찾을 것 입니다. 하지만 아직 PP-pathway 노드는 entrez ID 정수값 형태이므로 uniprot ID로 변경해 줄 필요가 있습니다. 이를 위해 편의상 파이썬 패키지로 ID 변환을 제공해주는 [PyEntrez]([https://pypi.org/project/pyEntrezId/](https://pypi.org/project/pyEntrezId/))를 사용합니다. 
  
`entrez2uniprot.py` 를 통해 entrez ID로 구성된 uniprot ID로 변환하여 파일로 저장할 수 있습니다. 이를 위해 PP-pathway gene id를 정리한 파일 `PP-Pathways_ppi.geneids.txt` 가 사용되며, 실행 결과 `PP-Pathways_ppi.entrez2uniprot.txt` 파일이 생성됩니다. 

    # PP-Pathways_ppi.entrez2uniprot.txt
    116685	P # Not found
    54069	Q9NYP9 # (entrez ID, uniprot ID)
    ...
    80208	NULL # Not found
    65201	Q9JK25 # (entrez ID, uniprot ID)
PP-pathway 데이터셋의 노든 노드에 대하여 변환과정을 수행한 결과, 대부분의 노드 ID는 변환되었으나 몇몇 노드는 그렇지 않았습니다. 우리는 변환된 노드 정보만 가지고 이후 데이터셋을 구성합니다
 

## 서열 기반 학습 및 평가 데이터셋 구성

PP-pathways edge list의 ID값을 UniProt ID 값에 해당하는 서열로 변환 후 학습용과 평가용으로 나누어 저장하는 작업을 `pp_pathways_sequences.py` 을 통해 얻을 수 있습니다. 
### (Input 1) UniProt fasta 포맷 서열 데이터
앞서 소개한 링크를 통해 다운 받을 수 있습니다

    # uniprot_sprot.fasta
    >sp|Q6GZX4|001R_FRG3G Putative transcription factor 001R OS=Frog virus 3 (isolate Goorha) OX=654924 GN=FV3-001R PE=4 SV=1
    MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPS
	EKGLIVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRRYQGIQTCKIPGKVLSDLD
	AKIKAYNLTVEGVEGFVRYSRVTKQHVAAFLKELRHSKQYENVNLIHYILTDKRVDIQHL
	EKDLVKDFKALVESAHRMRQGHMINVKYILYQLLKKHGHGPDGPDILTVKTGSKGVLYDD
	SFRKIYTDLGWKFTPL
	>sp|Q6GZX3|002L_FRG3G Uncharacterized protein 002L OS=Frog virus 3 (isolate Goorha) 	OX=654924 GN=FV3-002L PE=4 SV=1
	MSIIGATRLQNDKSDTYSAGPCYAGGCSAFTPRGTCGKDWDLGEQTCASGFCTSQPLCAR
	IKKTQVCGLRYSSKGKDPLVSAEWDSRGAPYVRCTYDADLIDTQAQVDQFVSMFGESPSLenter code here
	...

### (Input 2) Entrez ID 와 UniProt ID 쌍 
    # PP-Pathways_ppi.entrez2uniprot.txt
    116685	P # Not found
    54069	Q9NYP9 # (entrez ID, uniprot ID)
    ...
    80208	NULL # Not found
    65201	Q9JK25 # (entrez ID, uniprot ID)

### (Input 3) Human PPI edge list

    # PP-pathways_ppi.csv
    1394,2778 # pair-1(gene ID, gene ID)
    6331,17999 # pair-2(gene ID, gene ID)
    ...
    4790,79155 # pair-N(gene ID, gene ID)

### (Output) 학습 및 평가용으로 분리된 서열 쌍 데이터 셋
위 3가지 입력 데이터를 통해 네트워크 데이터셋에서 ID를 변환하고 각 ID에 대응하는 서열을 찾아 매칭시킵니다. 그리고 현재 edge list의 수가 너무 많기 때문에 약 30000개만 사용하며 ID값을 찾을 수 없는 일부 edge 정보들은 무시합니다. 그리고 서열 길이 역시 1000 이하만 사용하는 것으로 한다. 그 결과 학습 데이터셋 `PP-Pathways_ppi.train.txt`에는 17560개의 서열 쌍 그리고 각 노드별 서열 파일 `PP-Pathways_ppi.train.txt.seq` 에는 8503개의 서로다른 단백질 서열들이 존재하며, 평가용 데이터셋 `PP-Pathways_ppi.test.txt`에는 4390개 서열 쌍과 각 노드별 서열 파일 `PP-Pathways_ppi.test.txt.seq`에는 4366개 서로다른 단백질 서열들이 저장됩니다. 참고로 이 결과는 random seed에 따라 달라지기 때문에 대략적인 수치의 참고용 입니다.

## Position Specific Score Matrix(PSSM)
### PSI-BLAST 설치
단백질 서열을 사용한 대부분의 연구에 필수적으로 사용되는 특징 값 중 하나로 PSSM이 있습니다. PSSM은 쉽게 말해 단백질 서열의 profile 정보 입니다. 일반적으로 BLAST(Basic Local Alignment Search Tool)와 사용자가 지정한 서열 데이터베이스 파일에서 얻을 수 있습니다. 여기서는 PSI-BLAST(Position-Specific Iterated BLAST)를 Ubuntu 환경에서 설치하여 사용하였으며, 서열 데이터셋은 PDB를 사용합니다. 찾아야하는 PSSM 정보가 많기 때문에 웹 서비스 PSI-BLAST가 아닌 로컬 환경에 설치하여 사용해야 합니다.  PSI-BLAST 및 각종 BLAST 툴은 아래 명령어를 통해 쉽게 설치 할 수 있습니다

    sudo apt-get install ncbi-blast+


### BLAST 단백질(아미노산) 서열 데이터베이스 만들기
BLAST 설치가 완료되었다면 다음 할 일은 BLAST가 참조할 수 있는 데이터베이스를 구축하는 것 입니다. 데이터베이스 구축을 위해서는 단백질 서열 뭉치 데이터셋이 필요합니다. 데이터 베이스 구축 방법에 대한 내용은 [NCBI-BLAST 데이터베이스 도움말]([https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html](https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html))을 참고하시기 바랍니다. 해당 프로젝트에서는 [NCBI-BLAST FTP 서버]([https://ftp.ncbi.nlm.nih.gov/blast/db/](https://ftp.ncbi.nlm.nih.gov/blast/db/))에서 Protein Data Bank(PDB) [서열 뭉치](https://ftp.ncbi.nlm.nih.gov/blast/db/pdbaa.tar.gz)를 받아 데이터베이스로 변환 후 사용합니다. PDB 데이터베이스는 3차원 구조가 알려진 단백질들의 서열만 포함되어 있어 다른 데이터베이스보다 그 양이 작습니다. 우선 빠르게 PSSM을 찾기위해 PDB 데이터베이스를 사용합니다. 이후 스크립트는 `/ncbi/db/` 폴더에 변환된 PDB 데이터베이스를 저장하였음을 가정하고 진행합니다.

### 훈련 및 평가용 데이터셋의 노드 서열의 PSSM 구하기
이제 PSSM을 구하기 위한 설정은 모두 완료되었습니다. 편의를 위해 `seq2pssm.py` 는 멀티프로세싱을 통해 PSI-BLAST를 동작시키고 PSSM을 저장하는 작업을 수행할 수 있습니다.  별도의 입력 파라미터 설정은 없으므로 코드 수정 후 실행 바랍니다. 아래는 PSSM을 구하기 위한 PSI-BLAST 커맨드 예시 입니다.

    psiblast 
    -query input_sequence.fasta # 단일 서열 fasta 파일 
    -db /ncbi/db/pdbaa # PSSM 계산을 위한 PDB 데이터베이스 위치 
    -num_iterations 10 # 반복 회수 (자세한 내용은 도움말 참고) 
    -out_pssm out_pssm.txt # PSSM 연산 결과 저장용
    -out_ascii_pssm out.pssm # 이후 사용하게 될 (20 x N) 크기 PSSM이 저장되는 파일 
    -out out.txt # 출력 로그 생략용
   
   앞서 소개한 `seq2pssm.py` 를 사용하여 PSSM을 계산 할 수 있습니다.  단 PDB 데이터베이스에 포함되지 않는 서열들은 PSSM을 얻을 수 없기 때문에 이후 데이터 처리 과정에서 PSSM이 없는 서열 데이터는 생략합니다.
  
## 최종 데이터셋 구성
이제 서열 데이터 획득 및 학습 및 평가용 분리 그리고 PSSM 계산까지 완료하였습니다. 남은 작업은 해당 데이터를 통합하여 쉽게 사용할 수 있도록 단일 데이터베이스화 시키는 것 입니다. 이를 위해 numpy 패키지로 압축하여 저장합니다. 통합 과정을 `pssm2feature.py` 으로 수행 할 수 있으며, 결과 파일의 용량이 크기때문에 미리 하드디스크 공간을 넉넉히 비운 뒤 실행바랍니다. 

 

### 참고자료
1. Large-scale analysis of disease pathways in the human interactome. Monica Agrawal, Marinka Zitnik, and Jure Leskovec. Pacific Symposium on Biocomputing. 2018.
2. https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html
3. https://github.com/lwgray/pyEntrezId