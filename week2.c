#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define BUFSIZE 1024 //ファイルから読み込む一行の最大文字数
#define MAX_SEQ_NUM 30 //一つの転写因子に対して与えられる結合部位配列の最大数
#define MAX_GENE_NUM 8 /*与えられるプロモータ領域の最大遺伝子数*/

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

int main(int argc, char* argv[]){
  int seq_num = read_multi_seq(argv[1]); //１番目の引数で指定した転写因子の複数の結合部位配列を読み込む

  //頻度表の行列を定義する
  int hindo[4][BUFSIZE] = {0};  

//長さの検出
int max_len =0;
for(int i=0;i<BUFSIZE;i++){
  if(g_motif[0][i] != 0){
   max_len = max_len +1;
  }
  else{
    break;
  }
}


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
printf("hindhyo\n");
printf("A ");
for (int j = 0; j < max_len; j++) {
    printf("%3d ", hindo[0][j]);
    }
    printf("\n");

printf("C ");
for (int j = 0; j < max_len; j++) {
    printf("%3d ", hindo[1][j]);
    }
    printf("\n");

printf("G ");
for (int j = 0; j < max_len; j++) {
    printf("%3d ", hindo[2][j]);
    }
    printf("\n");

printf("T ");
for (int j = 0; j < max_len; j++) {
    printf("%3d ", hindo[3][j]);
    }
    printf("\n");

//対数オッズスコア行列の計算
int M ;
M = hindo[1][1]+hindo[2][1]+hindo[3][1]+hindo[4][1]+4;
float p[4][BUFSIZE] = {0};  
for(int i=0 ; i<4 ; i++){
  for(int j = 0; j < max_len; j++){
    p[i][j] = ((float)(hindo[i][j]+1)/(float)M);
  }
}

//手順1の結果の表示
printf("kakuritsu\n");
printf("A ");
for (int j = 0; j < max_len; j++) {
    printf("%f  ", p[0][j]);
    }
    printf("\n");

printf("C ");
for (int j = 0; j < max_len; j++) {
    printf("%f  ", p[1][j]);
    }
    printf("\n");

printf("G ");
for (int j = 0; j < max_len; j++) {
    printf("%f  ", p[2][j]);
    }
    printf("\n");

printf("T ");
for (int j = 0; j < max_len; j++) {
    printf("%f  ", p[3][j]);
    }
    printf("\n");

//手順2
float qx[4][BUFSIZE] = {0};
for(int i = 0; i<BUFSIZE ; i++){
qx[0][i] = (float)7519429/((float)7519429+(float)4637676+(float)7519429+(float)4637676);
qx[1][i] = (float)4637676/((float)7519429+(float)4637676+(float)7519429+(float)4637676);
qx[2][i] = (float)4637676/((float)7519429+(float)4637676+(float)7519429+(float)4637676);
qx[3][i] = (float)7519429/((float)7519429+(float)4637676+(float)7519429+(float)4637676);
}

printf("qx\n");
printf("A ");
for (int j = 0; j < max_len; j++) {
    printf("%f  ", qx[0][j]);
    }
    printf("\n");

printf("C ");
for (int j = 0; j < max_len; j++) {
    printf("%f  ", qx[1][j]);
    }
    printf("\n");

printf("G ");
for (int j = 0; j < max_len; j++) {
    printf("%f  ", qx[2][j]);
    }
    printf("\n");

printf("T ");
for (int j = 0; j < max_len; j++) {
    printf("%f  ", qx[3][j]);
    }
    printf("\n");





//手順3
float si[4][BUFSIZE];
for(int i=0 ; i<4 ; i++){
  for(int j = 0; j < max_len; j++){
    si[i][j] = log((p[i][j])/(qx[i][j]));
  }
}

printf("odds score\n");

printf("A ");
for (int j = 0; j < max_len; j++) {
    printf("%10f  ", si[0][j]);
    }
    printf("\n");

printf("C ");
for (int j = 0; j < max_len; j++) {
    printf("%10f  ", si[1][j]);
    }
    printf("\n");

printf("G ");
for (int j = 0; j < max_len; j++) {
    printf("%10f  ", si[2][j]);
    }
    printf("\n");

printf("T ");
for (int j = 0; j < max_len; j++) {
    printf("%10f  ", si[3][j]);
    }
    printf("\n\n\n");



  int gene_num = read_promoter(argv[2]);  //２番目の引数で指定した遺伝子のプロモータ領域を読み込む
  
  printf("promoter_sequence:\n");
  for(int i = 0; i < gene_num; i++){
    printf(">%s\n", g_pro[i].name); //読み込んだプロモータ領域を表示
    printf("%s\n", g_pro[i].seq);
  }

int len_seq =0;
for(int i=0;i<BUFSIZE;i++){
  if(g_pro[0].seq[i] != '\0'){
   len_seq = len_seq +1;
  }
  else{
    break;
  }
}

int hit,max_hit;
hit=0;
max_hit=0;







printf("%d\n",len_seq);









  return 0;
}

