#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define BUFSIZE 1024 //ファイルから読み込む一行の最大文字数
#define MAX_SEQ_NUM 30 //一つの転写因子に対して与えられる結合部位配列の最大数
#define MAX_GENE_NUM 8 /*与えられるプロモータ領域の最大遺伝子数*/
#define BASES_NUM 4
#define threshold 4 //閾値

char g_motif[MAX_SEQ_NUM][BUFSIZE]; //転写因子の結合部位配列を保存する配列

struct promoter{
  char name[BUFSIZE];
  char seq[BUFSIZE];
}g_pro[MAX_GENE_NUM]; //遺伝子のプロモータ領域を保存する構造体

//グローバル変数はローカル変数と区別するため、全部大文字にするかg_を先頭につけるのが一般的

int read_multi_seq(char* filename){
  int seq_num = 0;
  char buffer[BUFSIZE];
  FILE *fp = fopen(filename,"r");

  if(fp == NULL){
    printf("motif_region_file open error.\n");
    exit(1); //ファイルが開けなかった場合プログラムを終了
  }

  while(fscanf(fp, "%s", buffer) != EOF){ //プログラムから一行ずつ読み込む
    if(buffer[strlen(buffer)-1]=='\n'){
      buffer[strlen(buffer)-1]='\0'; //改行を切り落とす
    }
    strcpy(g_motif[seq_num],buffer); //結合部位配列を保存
    seq_num++;
  }
  return seq_num;
}

int read_promoter(char *filename){
  int gene_num = 0;  
  char buffer[BUFSIZE];
  FILE *fp = fopen(filename,"r");

  if(fp == NULL){
    printf("scorefile open error.\n");
    exit(1);
  }

  while(fscanf(fp, "%s", buffer) != EOF){
    if(buffer[strlen(buffer)-1]=='\n'){
      buffer[strlen(buffer)-1]='\0';
    }
    
    if(buffer[0]=='>'){
      strcpy(g_pro[gene_num].name,buffer+1); 
    }else{
      strcpy(g_pro[gene_num].seq,buffer);
      gene_num++;
    }    
  }
  return gene_num;
}

int get_max_len(char *motif){
  int len = 0;
  while (motif[len] != '\0') {
        len++;
    }
  return len;
}

void print_matrix_int(int matrix[4][BUFSIZE], int max_len) {
    char enki[4] = {'A', 'C', 'G', 'T'};
    for (int i = 0; i < 4; i++) {
        printf("%c ", enki[i]);

        for (int j = 0; j < max_len; j++) {
            printf("%3d", matrix[i][j]);
        }
        printf("\n");
    }
}


void print_matrix(float matrix[4][BUFSIZE], int max_len) {
    char enki[BASES_NUM] = {'A', 'C', 'G', 'T'};
    for (int i = 0; i < BASES_NUM; i++) {
        printf("%c ", enki[i]);

        for (int j = 0; j < max_len; j++) {
            printf("%10f", matrix[i][j]);
        }
        printf("\n");
    }
}



//以下メイン関数
int main(int argc, char* argv[]){
  int seq_num = read_multi_seq(argv[1]); //１番目の引数で指定した転写因子の複数の結合部位配列を読み込む

//頻度表の行列を定義する
int hindo[4][BUFSIZE] = {0};  

//長さの検出
int max_len = get_max_len(g_motif[0]);

// 頻度表を作る
for (int col = 0; col < max_len; col++) {
    for (int row = 0; row < seq_num; row++) {
        char c = g_motif[row][col];
        switch (c) {
            case 'A':
                hindo[0][col]++;
                break;
            case 'C':
                hindo[1][col]++;
                break;
            case 'G':
                hindo[2][col]++;
                break;
            case 'T':
                hindo[3][col]++;
                break;
        }
    }
}

//頻度表の表示
//printf("hindhyo\n");
//print_matrix_int(hindo,max_len);

//手順1の計算
int M ;
M = hindo[0][1]+hindo[1][1]+hindo[2][1]+hindo[3][1]+4;
float p[BASES_NUM][BUFSIZE] = {0};  
for(int i=0 ; i<BASES_NUM ; i++){
  for(int j = 0; j < max_len; j++){
    p[i][j] = ((float)(hindo[i][j]+1)/(float)M);
  }
}

//手順1の結果の表示
//printf("kakuritsu\n");
//print_matrix(p, max_len);
//printf("\n");

//手順2の計算
float qx[BASES_NUM][BUFSIZE] = {0};
for(int i = 0; i<BUFSIZE ; i++){
qx[0][i] = (float)7519429/((float)7519429+(float)4637676+(float)7519429+(float)4637676);
qx[1][i] = (float)4637676/((float)7519429+(float)4637676+(float)7519429+(float)4637676);
qx[2][i] = (float)4637676/((float)7519429+(float)4637676+(float)7519429+(float)4637676);
qx[3][i] = (float)7519429/((float)7519429+(float)4637676+(float)7519429+(float)4637676);
}

//手順2の表示
//printf("qx\n");
//print_matrix(qx, max_len);
//printf("\n");


//手順3の計算
float si[BASES_NUM][BUFSIZE];
for(int i=0 ; i<4 ; i++){
  for(int j = 0; j < max_len; j++){
    si[i][j] = log((p[i][j])/(qx[i][j]));
  }
}

//手順3の表示
//printf("odds score\n");
//print_matrix(si, max_len);
//printf("\n\n\n");


//
//以下プロモーター配列との比較
//

int gene_num = read_promoter(argv[2]);  //２番目の引数で指定した遺伝子のプロモータ領域を読み込む

//seqの長さを確認しlen_seqへ代入
int len_seq = get_max_len(g_pro[0].seq);

//必要な変数の定義
float hit=0;
int overnum[BUFSIZE] = {0};
int number;
float overhit[BUFSIZE] = {0};

printf("%s\n\n", argv[1]);

//プロモーターを回すfor文
for(int j = 0;j<gene_num;j++){
  printf("pro:%s\n\n", g_pro[j].name);
number = 0;

//それぞれのプロモーター配列におけるスコアを算出
for(int t=0;t<len_seq-max_len;t++){
    hit=0;
    for(int i=t;i<t+max_len;i++){
        char c = g_pro[j].seq[i];
        switch (c) {
            case 'A':
                hit = hit + si[0][i-t];
                break;
            case 'C':
                hit = hit + si[1][i-t];
                break;
            case 'G':
                hit = hit + si[2][i-t];
                break;
            case 'T':
                hit = hit + si[3][i-t];
                break;
        }
}

//閾値の判定，記録
if(hit > threshold){
            overhit[number]= hit;
            overnum[number]= t;
            number = number +1;
        }
}

//諸々の表示
if (number == 0) {
        printf("none\n\n");
    } else {
for(int b = 0; b < number; b++) {
   printf("pos:%d",overnum[b]);
   printf("\n");
   for(int l = overnum[b]; l < overnum[b] + max_len; l++) {
        printf("%c", g_pro[j].seq[l]);
    }
   printf("\n");
   printf("hit : %f",overhit[b]);
   printf("\n\n");

}
printf("\n");
}

}
return 0;

}