#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define BUFSIZE 1024 //ファイルから読み込む一行の最大文字数
#define MAX_SEQ_NUM 30 //一つの転写因子に対して与えられる結合部位配列の最大数
#define MAX_GENE_NUM 8 /*与えられるプロモータ領域の最大遺伝子数*/
#define BASES_NUM 4
#define threshold 4 //閾値
#define SIZE 500
#define NUM 50

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


int get_max_len(char *motif){
  int len = 0;
  while (motif[len] != '\0') {
        len++;
    }
  return len;
}

//ランダム配列の取得
void get_genome(char genome[NUM][SIZE]){
     int seq[SIZE];
    srand((unsigned) time(NULL));
     for (int i = 0; i < SIZE; i++) {
        seq[i] = rand() % 10000; // 0〜9999の範囲で乱数
    }

    for(int t=0;t<NUM;t++){
    for(int i=0;i<SIZE;i++){
        if (0 <= seq[i] && seq[i] <= 3092){
            genome[t][i]='A';
        }
        else if(3093 <= seq[i] && seq[i] <= 5000){
            genome[t][i]='C';
        }
        else if((5001 <= seq[i] && seq[i] <= 6908)){
            genome[t][i]='G';
        }
        else{
            genome[t][i]='T';
        }
    }
}
}

//ランダム配列の表示
void matrix(char genome1[NUM][SIZE]){
    for(int t=0;t<NUM;t++){
        for(int i=0;i<SIZE;i++){
            printf("%c",genome1[t][i]);
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


//手順1の計算
int M ;
M = hindo[0][1]+hindo[1][1]+hindo[2][1]+hindo[3][1]+4;
float p[BASES_NUM][BUFSIZE] = {0};  
for(int i=0 ; i<BASES_NUM ; i++){
  for(int j = 0; j < max_len; j++){
    p[i][j] = ((float)(hindo[i][j]+1)/(float)M);
  }
}



//手順2の計算
float qx[BASES_NUM][BUFSIZE] = {0};
for(int i = 0; i<BUFSIZE ; i++){
qx[0][i] = (float)7519429/((float)7519429+(float)4637676+(float)7519429+(float)4637676);
qx[1][i] = (float)4637676/((float)7519429+(float)4637676+(float)7519429+(float)4637676);
qx[2][i] = (float)4637676/((float)7519429+(float)4637676+(float)7519429+(float)4637676);
qx[3][i] = (float)7519429/((float)7519429+(float)4637676+(float)7519429+(float)4637676);
}




//手順3の計算
float si[BASES_NUM][BUFSIZE];
for(int i=0 ; i<4 ; i++){
  for(int j = 0; j < max_len; j++){
    si[i][j] = log((p[i][j])/(qx[i][j]));
  }
}


char genome[NUM][SIZE];
    get_genome(genome);

float score[NUM][SIZE] = {0};

//seqの長さを確認しlen_seqへ代入
int len_seq = get_max_len(g_pro[0].seq);

//必要な変数の定義
float hit=0;
int overnum[BUFSIZE] = {0};
int number;
float overhit[BUFSIZE] = {0};

//それぞれのプロモーター配列におけるスコアを算出
for(int t=0; t<NUM; t++){
    for(int i=0; i<SIZE - max_len; i++){
        float hit = 0;
        for(int j=0; j<max_len; j++){
            char c = genome[t][i+j];
            switch (c) {
                case 'A': hit += si[0][j]; 
                break;
                case 'C': hit += si[1][j]; 
                break;
                case 'G': hit += si[2][j]; 
                break;
                case 'T': hit += si[3][j]; 
                break;
            }
        }
        score[t][i] = hit;
    }
}


FILE *fp = fopen("scorefile.txt", "w");

for(int i=0;i<NUM;i++){
    for(int j=0;j<SIZE;j++){
        fprintf(fp,"%f",score[i][j]);
    }
    fprintf(fp,"\n\n");
}

return 0;
}