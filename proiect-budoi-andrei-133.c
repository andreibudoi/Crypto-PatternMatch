#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h> // Librarie pentru tipurile de int (uint8_t, uint16_t, ect..)
#pragma pack(1) // Formatare Struct Header (Sterge paddingul struct-ului
// pentru a citi corect informatiile din header)

//VALORILE PENTRU INCADRAREA FIECAREI CIFRE
#define CULOARE_0 255,0,0
#define CULOARE_1 255,255,0
#define CULOARE_2 0,255,0
#define CULOARE_3 0,255,255
#define CULOARE_4 255,0,255
#define CULOARE_5 0,0,255
#define CULOARE_6 192,192,192
#define CULOARE_7 255,140,0
#define CULOARE_8 128,0,128
#define CULOARE_9 128,0,0

//Structura pixel ce retine valorile pentru fiecare culoare pentru un .bmp de tip 24bit
typedef struct pixelBGR
{
    //fiecare valoare este retinuta pe un octet
    uint8_t blue;
    uint8_t green;
    uint8_t red;

} pixelBGR;

//Structura ce retine informatiile din header
typedef struct HeaderInfo
{
    //Detaliile despre header sunt gasite pe Wikipedia
    uint16_t Field;
    uint32_t BMPSize;
    uint16_t Reserved1;
    uint16_t Reserved2;
    uint32_t Offset;
    //////////////////
    uint32_t InfoSize;
    int32_t Width;
    int32_t Height;
    uint16_t ColorPlanes;
    uint16_t ColorDepth;
    uint32_t Compression;
    uint32_t ImageSize;
    int32_t VerticalRez;
    int32_t HorizontalRez;
    uint32_t Colors;
    uint32_t Important;

} HeaderInfo;

typedef struct detectie{
    uint32_t x,y;
    double corr;
    uint8_t red,green,blue;
}detectie;

//Functie ce genereaza un vector cu numere pseudo-random folosind algoritmul XorShift32
uint32_t *XorShift(uint32_t seed,uint32_t W, uint32_t H)
{
    uint32_t i,r=seed;

    uint32_t *xorarray=(uint32_t*)malloc(2*W*H*sizeof(uint32_t));
    xorarray[0]=seed;
    for(i=1; i<2*W*H; i++)
    {
        r ^= r << 13;
        r ^= r >> 17;
        r ^= r << 5;
        xorarray[i]=r;
    }
    return xorarray;
}

//Functie pentru testul chi patrat
void chi (HeaderInfo header, pixelBGR *forma_liniara,char *numepoza)
{
    double chiR=0,chiG=0,chiB=0,f_barat;
    int *frecR,*frecG,*frecB;
    frecR=(int *)calloc(256,sizeof(int));
    frecG=(int *)calloc(256,sizeof(int));
    frecB=(int *)calloc(256,sizeof(int));
    f_barat=(header.Width)*(header.Height)/256;
    int i;
    for(i=0; i<(header.Width)*(header.Height); i++)
    {
        frecR[forma_liniara[i].red]++;
        frecG[forma_liniara[i].green]++;
        frecB[forma_liniara[i].blue]++;

    }
    for(i=0; i<256; i++)
    {
        chiR=chiR+((frecR[i]-f_barat)*(frecR[i]-f_barat))/f_barat;
        chiG=chiG+((frecG[i]-f_barat)*(frecG[i]-f_barat))/f_barat;
        chiB=chiB+((frecB[i]-f_barat)*(frecB[i]-f_barat))/f_barat;
    }

    fprintf(stdout,"Chi-squared test on RGB channels for %s:\n",numepoza);

    fprintf(stdout,"R: ");
    fprintf(stdout,"%.2f",chiR);
    fprintf(stdout,"\n");

    fprintf(stdout,"G: ");
    fprintf(stdout,"%.2f",chiG);
    fprintf(stdout,"\n");

    fprintf(stdout,"B: ");
    fprintf(stdout,"%.2f",chiB);
    fprintf(stdout,"\n");
    fprintf(stdout,"\n");


    free(frecR);
    free(frecG);
    free(frecB);
}

void greyscale(pixelBGR **matrice,HeaderInfo header)
{
    uint32_t i,j;
    uint8_t aux=0;

    for(i=0;i<header.Height;i++)
        for(j=0;j<header.Width;j++)
        {
            aux=0.299* matrice[i][j].red + 0.587 * matrice[i][j].green + 0.114 * matrice[i][j].blue;
            matrice[i][j].red=matrice[i][j].green=matrice[i][j].blue=aux;
        }
}

//Citeste R0 si SV
void f_seed(uint32_t *xorseed, uint32_t *SV, char secretKey[])
{
    FILE *f;
    f = fopen(secretKey, "r");
    fscanf(f, "%d %d", &(*xorseed), &(*SV));
    fclose(f);
}

//Verifica daca fisierul a fost deschis corect
void FileErrorCheck(FILE *f,char *numepoza)
{
    if (f == NULL)
    {
        fprintf(stderr, "Fisierul %s nu a putut fi deschis :( ", numepoza);
        exit(0);
    }
    printf("Fisierul %s a fost deschis cu succes\n", numepoza);
}

//functie ce modifica culoarea pixelilor
pixelBGR *modificarePixeli(HeaderInfo header, pixelBGR *forma_liniara,uint32_t SV,uint32_t xorseed)
{
    uint32_t i, *xorarray;
    xorarray=XorShift(xorseed,header.Width,header.Height);

    pixelBGR *forma_liniara_enc=(pixelBGR *)malloc((header.Width)*(header.Height)*sizeof(pixelBGR));

    //Modificarea primului pixel:

    //BLUE
    forma_liniara_enc[0].blue=forma_liniara[0].blue ^ (SV&255)^ (xorarray[(header.Width)*(header.Height)]&255);

    //shiftare spre dreapta cu un byte
    SV=SV>>8;
    xorarray[(header.Width)*(header.Height)]=xorarray[(header.Width)*(header.Height)]>>8;

    //GREEN
    forma_liniara_enc[0].green=(forma_liniara[0].green) ^ (SV&255)^ (xorarray[(header.Width)*(header.Height)]&255);

    //shiftare spre dreapta cu un byte
    SV=SV>>8;
    xorarray[(header.Width)*(header.Height)]=xorarray[(header.Width)*(header.Height)]>>8;

    //RED
    forma_liniara_enc[0].red=(forma_liniara[0].red) ^ (SV&255)^ (xorarray[(header.Width)*(header.Height)]&255);

    //*/////////////*//

    //Modificarea urmatorilor k pixeli:
    for(i=1; i<(header.Width)*(header.Height); i++)
    {
        //BLUE
        forma_liniara_enc[i].blue=forma_liniara[i].blue ^ forma_liniara_enc[i-1].blue ^ (xorarray[(header.Width)*(header.Height)+i]&255);

        //shiftare spre dreapta cu un byte
        xorarray[(header.Width)*(header.Height)+i]=xorarray[(header.Width)*(header.Height)+i]>>8;

        //GREEN
        forma_liniara_enc[i].green=(forma_liniara[i].green) ^ (forma_liniara_enc[i-1].green) ^ (xorarray[(header.Width)*(header.Height)+i]&255);

        //shiftare spre dreapta cu un byte
        xorarray[(header.Width)*(header.Height)+i]=xorarray[(header.Width)*(header.Height)+i]>>8;

        //RED
        forma_liniara_enc[i].red=(forma_liniara[i].red) ^ (forma_liniara_enc[i-1].red) ^ (xorarray[(header.Width)*(header.Height)+i]&255);
    }

    free(xorarray);
    return forma_liniara_enc;
}

pixelBGR *modificarePixeli_decript(HeaderInfo header, pixelBGR *forma_liniara_enc,uint32_t SV,uint32_t xorseed)
{
    uint32_t i, *xorarray;
    xorarray=XorShift(xorseed,header.Width,header.Height);

    pixelBGR *forma_liniara=(pixelBGR *)malloc((header.Width)*(header.Height)*sizeof(pixelBGR));

    //Modificarea primului pixel:

    //BLUE
    forma_liniara[0].blue=forma_liniara_enc[0].blue ^ (SV&255)^ (xorarray[(header.Width)*(header.Height)]&255);

    //Modificarea primului pixel:
    SV=SV>>8;
    xorarray[(header.Width)*(header.Height)]=xorarray[(header.Width)*(header.Height)]>>8;

    //GREEN
    forma_liniara[0].green=(forma_liniara_enc[0].green) ^ (SV&255)^ (xorarray[(header.Width)*(header.Height)]&255);

    //Modificarea primului pixel:
    SV=SV>>8;
    xorarray[(header.Width)*(header.Height)]=xorarray[(header.Width)*(header.Height)]>>8;

    //RED
    forma_liniara[0].red=(forma_liniara_enc[0].red) ^ (SV&255)^ (xorarray[(header.Width)*(header.Height)]&255);

    //*/////////////*//

    //Modificarea urmatorilor k pixeli:
    for(i=1; i<(header.Width)*(header.Height); i++)
    {
        //BLUE
        forma_liniara[i].blue=forma_liniara_enc[i].blue ^ forma_liniara_enc[i-1].blue ^ (xorarray[(header.Width)*(header.Height)+i]&255);

        //Modificarea primului pixel:
        xorarray[(header.Width)*(header.Height)+i]=xorarray[(header.Width)*(header.Height)+i]>>8;

        //GREEN
        forma_liniara[i].green=(forma_liniara_enc[i].green) ^ (forma_liniara_enc[i-1].green) ^ (xorarray[(header.Width)*(header.Height)+i]&255);

        //Modificarea primului pixel:
        xorarray[(header.Width)*(header.Height)+i]=xorarray[(header.Width)*(header.Height)+i]>>8;

        //RED
        forma_liniara[i].red=(forma_liniara_enc[i].red) ^ (forma_liniara_enc[i-1].red) ^ (xorarray[(header.Width)*(header.Height)+i]&255);
    }

    free(xorarray);
    return forma_liniara;
}

//Functie ce permuta pixelii din matricea liniariata
pixelBGR *permutare(HeaderInfo header, pixelBGR *forma_liniara, uint32_t xorseed)
{
    uint32_t i, j, *xorarray, *P,aux;
    xorarray=XorShift(xorseed,header.Width,header.Height);
    P=(uint32_t*)malloc((header.Width)*(header.Height)*sizeof(uint32_t));

    //Permutare initiala
    for(i=0; i<(header.Width)*(header.Height); i++)
        P[i]=i;

    pixelBGR *forma_liniara_perm=(pixelBGR *)malloc((header.Width)*(header.Height)*sizeof(pixelBGR));

    //Creare permutare randomizata folosing algoritmul lui Durstenfeld si XorShift32
    int k=1;
    for (i = (header.Width)*(header.Height) - 1; i >= 1; i--)
    {
        j = xorarray[k++] % (i+1);

        aux=P[j];
        P[j]=P[i];
        P[i]=aux;
    }


    //Permutarea pixelilor
    for (i =0; i <(header.Width)*(header.Height); i++)
        forma_liniara_perm[P[i]]=forma_liniara[i];



    free(xorarray);
    free(P);
    return forma_liniara_perm;


}
pixelBGR *permutare_decript(HeaderInfo header, pixelBGR *forma_liniara_perm, uint32_t xorseed)
{
    uint32_t i, j, *xorarray, *P,aux;
    xorarray=XorShift(xorseed,header.Width,header.Height);
    P=(uint32_t*)malloc((header.Width)*(header.Height)*sizeof(uint32_t));

    //Permutare initiala
    for(i=0; i<(header.Width)*(header.Height); i++)
        P[i]=i;

    pixelBGR *forma_liniara=(pixelBGR *)malloc((header.Width)*(header.Height)*sizeof(pixelBGR));

    //Creare permutare randomizata folosing algoritmul lui Durstenfeld si XorShift32
    int k=1;
    for (i = (header.Width)*(header.Height) - 1; i >= 1; i--)
    {
        j = xorarray[k++] % (i+1);

        aux=P[j];
        P[j]=P[i];
        P[i]=aux;
    }

    //Permutarea pixelilor
    for (i =0; i <(header.Width)*(header.Height); i++)
        forma_liniara[i]=forma_liniara_perm[P[i]];



    free(xorarray);
    free(P);
    return forma_liniara;


}

pixelBGR *citire_encript(char numepoza[], HeaderInfo *header)
{
    FILE *f;
    f = fopen(numepoza, "rb"); //Deschide fisierul
    FileErrorCheck(f,numepoza); // Verificare fisier

    //Citire header
    fread(&(*header), sizeof(HeaderInfo), 1, f);

    //Alocare memorie
    pixelBGR *forma_liniara;
    forma_liniara = (pixelBGR*)malloc((header->Width) * (header->Height) * sizeof(pixelBGR));

    //Calculare padding
    int padding;
    if (header->Width % 4 != 0)
        padding = 4 - (3 * header->Width) % 4;
    else
        padding = 0;

    int i, j;
    //Liniarizarea pixelilor citind un array de pixeli
    for (i = 0; i < header->Height; i++)
    {
        for (j = 0; j < header->Width; j++)
        {
            fread(&forma_liniara[(header->Height - 1 - i) *  header->Width + j], sizeof(pixelBGR), 1, f);
        }
        //Sar peste padding
        fseek(f, padding, SEEK_CUR);
    }

    return forma_liniara;
    fclose(f);
}

void afisare_encript(char numepoza_mod[], HeaderInfo header, pixelBGR *forma_liniara)
{
    FILE *f;
    f = fopen(numepoza_mod, "wb"); //Deschide fisierul
    FileErrorCheck(f,numepoza_mod); //Verificare fisier
    char zero = '\0'; //Byte gol pentru padding

    //Scriere header
    fwrite(&header, sizeof(HeaderInfo), 1, f);

    //Calculare padding
    int padding;
    if (header.Width % 4 != 0)
        padding = 4 - (3 * header.Width) % 4;
    else
        padding = 0;

    int i, j;
    //Scriere pixeli in fisier
    for (i = 0; i < header.Height; i++)
    {
        for (j = 0; j < header.Width; j++)
        {
            fwrite(&forma_liniara[(header.Height-1-i) *  header.Width + j], sizeof(pixelBGR), 1, f);
        }
        //Se pune paddingul
        fwrite(&zero, padding, 1, f);
    }


    fclose(f);

}

void criptare(char numepoza[], char numepoza_mod[], char secretKey[])
{
    uint32_t xorseed, SV;
    pixelBGR *forma_liniara,*forma_liniara_perm,*forma_liniara_enc;
    HeaderInfo header;

    f_seed(&xorseed, &SV, secretKey);

    forma_liniara = citire_encript(numepoza, &header);		//citire poza de criptat
    forma_liniara_perm = permutare(header, forma_liniara, xorseed); //permutarea pixelilor
    forma_liniara_enc = modificarePixeli(header, forma_liniara_perm, SV, xorseed); //modificarea culorilor folosind codul secret

    //crearea pozei criptate
    afisare_encript(numepoza_mod, header, forma_liniara_enc);
    printf("Poza a fost criptata cu succes!\n");

    //testul chi patrat
    chi(header,forma_liniara,numepoza);
    chi(header,forma_liniara_enc,numepoza_mod);

    //eliberarea memoriei
    free(forma_liniara_enc);
    free(forma_liniara_perm);
    free(forma_liniara);
}

void decriptare(char numepoza[], char numepoza_mod[], char secretKey[])
{
    //declarari
    uint32_t xorseed, SV;
    pixelBGR *forma_liniara,*forma_liniara_perm,*forma_liniara_enc;
    HeaderInfo header;

    f_seed(&xorseed, &SV, secretKey); //citire cod secret

    forma_liniara_enc = citire_encript(numepoza_mod, &header); //citire poza de decriptat
    forma_liniara_perm = modificarePixeli_decript(header,forma_liniara_enc,SV,xorseed); //permutarea pixelilor
    forma_liniara = permutare_decript(header, forma_liniara_perm, xorseed); //modificarea culorilor folosind codul secret

    //crearea pozei decriptate
    afisare_encript(numepoza, header, forma_liniara);
    printf("Poza a fost decriptata cu succes!\n\n");

    //eliberarea memoriei
    free(forma_liniara_enc);
    free(forma_liniara_perm);
    free(forma_liniara);
}

/////////////////////////////////
/////////////////////////////////

//Functie ce calculeaza suprapunerea spatiala a doua ferestre
double intersect (uint32_t n, uint32_t L1,uint32_t R1, uint32_t B1 ,uint32_t U1,uint32_t L2,uint32_t R2, uint32_t B2,uint32_t U2)
{
    //left, right, bottom, upper
    int32_t L3,R3,B3,U3,dx,dy;

    if(L1<L2)L3=L2;
    else L3=L1;

    if(R1<R2) R3=R1;
    else R3=R2;

    if(B1<B2) B3=B1;
    else B3=B2;

    if(U1<U2) U3=U2;
    else U3=U1;
    dx=R3-L3;
    dy=B3-U3;
    int intersectArea;
    if (dx < 0 || dy <0) return 0;
    dx++;
    dy++;

    intersectArea=dx*dy;
    double val = (double)(intersectArea)/(double)(n+n-intersectArea);

    return val;
}

//Functie ce elimina detectiile inutile
void eliminare(detectie *D, uint32_t *k,HeaderInfo header)
{
    uint32_t i,j,l;
    double val;
    for(i=0;i<(*k)-1;i++)
        for(j=i+1;j<(*k);j++)
        {
            val=intersect(header.Height*header.Width,D[i].y - (header.Width/2), D[i].y + (header.Width/2),D[i].x +(header.Height/2),D[i].x -(header.Height/2),
                                                    D[j].y - (header.Width/2), D[j].y + (header.Width/2),D[j].x +(header.Height/2),D[j].x -(header.Height/2));

            if(val > 0.2)
            {

                for(l=j;l<(*k)-1;l++)
                    D[l]=D[l+1];
                j--;
                (*k)--;

            }
        }
}

//Functia de comparare pentru qsort
int cmp(const void* a, const void* b)
{
   detectie *A=(detectie*)a;
   detectie *B=(detectie*)b;
   if(B->corr > A->corr)return 1;
   else if(B->corr < A->corr) return -1;
   else return 0;
}

//Functie ce deseneaza conturul unei ferestre cu o culoare anume
void incadrare (pixelBGR **source_mod,HeaderInfo header_source,uint32_t H, uint32_t W,uint32_t x, uint32_t y, uint8_t R,  uint8_t G, uint8_t B)
{
    uint32_t i,j;

    for(i=( x - (H/2) ); i <= ( x + (H/2) )  ; i++)
    {
        source_mod[i][y - (W/2)].blue=B;
        source_mod[i][y - (W/2)].green=G;
        source_mod[i][y - (W/2)].red=R;

        source_mod[i][y + (W/2)].blue=B;
        source_mod[i][y + (W/2)].green=G;
        source_mod[i][y + (W/2)].red=R;
    }

    for(j=(y - (W/2)); j<= (y + (W/2)) ;j++)
    {
        source_mod[x - (H/2)][j].blue=B;
        source_mod[x - (H/2)][j].green=G;
        source_mod[x - (H/2)][j].red=R;

        source_mod[x + (H/2) ][j].blue=B;
        source_mod[x + (H/2) ][j].green=G;
        source_mod[x + (H/2) ][j].red=R;
    }
}

//Calculeaza corelatia dintre o fereastra si un sablon folosind algoritmul din cerinta
double corr(pixelBGR **source, pixelBGR **s ,HeaderInfo header_source,HeaderInfo header_sablon,uint32_t x, uint32_t y)
{

    double n=(header_sablon.Height*header_sablon.Width);
    double suma_corr=0;
    double medieSablon=0, medieFereastra=0;
    double sigmaSablon=0,sigmaFereastra=0;
    uint32_t i,j;

    //Se calculeaza media valorilor intensităților grayscale a pixelilor din sablon
    for(i=0 ; i< header_sablon.Height  ; i++)
        for(j=0; j< header_sablon.Width ;j++)
            medieSablon=medieSablon+s[i][j].red;
    medieSablon=medieSablon/(header_sablon.Height*header_sablon.Width);

    //Se calculeaza media valorilor intensităților grayscale a pixelilor din fereastra
    for(i=( x - (header_sablon.Height/2) )  ; i<= ( x + (header_sablon.Height/2) )  ; i++)
        for(j=(y - (header_sablon.Width/2)); j<= (y + (header_sablon.Width/2)) ;j++)
            medieFereastra=medieFereastra+source[i][j].red;
    medieFereastra=medieFereastra/n;

    //Se calculeaza deviația standard a valorilor intensităților grayscale a pixelilor din șablon
    for(i=0 ; i< header_sablon.Height; i++)
        for(j=0; j< header_sablon.Width ;j++)
            sigmaSablon=sigmaSablon+(s[i][j].red - medieSablon)*(s[i][j].red - medieSablon);
    sigmaSablon=sqrt(sigmaSablon * (1/(n-1)));

    //Se calculeaza deviația standard a valorilor intensităților grayscale a pixelilor din fereastra
    for(i=( x - (header_sablon.Height/2) )  ; i<= ( x + (header_sablon.Height/2) )  ; i++)
        for(j=(y - (header_sablon.Width/2)); j<= (y + (header_sablon.Width/2)) ;j++)
             sigmaFereastra=sigmaFereastra+(source[i][j].red - medieFereastra)*(source[i][j].red - medieFereastra);
    sigmaFereastra=sqrt(sigmaFereastra* (1/(n-1)));

    //Se calculeaza corelatia
    for(i=( x - (header_sablon.Height/2) )  ; i<= ( x + (header_sablon.Height/2) )  ; i++)
        for(j=(y - (header_sablon.Width/2)); j<= (y + (header_sablon.Width/2)) ;j++)
            suma_corr=suma_corr+( (1/(sigmaFereastra * sigmaSablon))* (source[i][j].red - medieFereastra)* (s[i - (x - header_sablon.Height/2)][j - (y - (header_sablon.Width/2))].red - medieSablon)  );
    suma_corr=suma_corr*(1/n);

    return suma_corr;
}

//Functie pentru template matching
void matching(pixelBGR **source, pixelBGR **s ,pixelBGR **source_mod,HeaderInfo header_source, HeaderInfo header_sablon, detectie *D,uint32_t *k,double pS, uint8_t R,  uint8_t G, uint8_t B)
{
    uint32_t i,j;
    float val_corr;

    // Se gliseaza sablonul s
    for(i=header_sablon.Height/2 ; i<( header_source.Height - (header_sablon.Height/2) ) ; i++)
        for(j=header_sablon.Width/2 ; j< ( header_source.Width - (header_sablon.Width/2)) ;j++)
            {
                //Se calculeaza corelatia
                val_corr=corr(source, s , header_source, header_sablon,i,j);
                if(val_corr>pS)
                {
                    //Se adauga in tabloul de detectii
                    D[*k].red=R;
                    D[*k].green=G;
                    D[*k].blue=B;
                    D[*k].corr=val_corr;
                    D[*k].x=i;
                    D[*k].y=j;
                    (*k)++;
                }
            }
}

//Se citeste imaginea intr-o matrice de pixeli
pixelBGR *citire_mat(char numepoza[], HeaderInfo *header)
{
    FILE *f;
    f = fopen(numepoza, "rb"); //Deschide fisierul
    FileErrorCheck(f,numepoza); // Verificare fisier

    //Citire header
    fread(&(*header), sizeof(HeaderInfo), 1, f);

    //Alocare memorie
    int i, j;
    pixelBGR **matrice_pixeli=(pixelBGR**)malloc((header->Height) * sizeof(pixelBGR*));
    for(i = 0; i < header->Height; i++) matrice_pixeli[i] = (int *)calloc(header->Width, sizeof(pixelBGR));

    //Calculare padding
    int padding;
    if (header->Width % 4 != 0)
        padding = 4 - (3 * header->Width) % 4;
    else
        padding = 0;


    //Liniarizarea pixelilor citind un array de pixeli
    for (i = 0; i < header->Height; i++)
    {
        for (j = 0; j < header->Width; j++)
        {
            fread(&matrice_pixeli[i][j], sizeof(pixelBGR), 1, f);
        }
        //Sar peste padding
        fseek(f, padding, SEEK_CUR);
    }

    return matrice_pixeli;
    fclose(f);
}

//Se creaza o noua poza folosind matricea de pixeli
void afisare_mat(char numepoza_mod[], HeaderInfo header, pixelBGR **matrice_pixeli)
{
    FILE *f;
    f = fopen(numepoza_mod, "wb"); //Deschide fisierul
    FileErrorCheck(f,numepoza_mod); //Verificare fisier
    char zero = '\0'; //Byte gol pentru padding

    //Scriere header
    fwrite(&header, sizeof(HeaderInfo), 1, f);

    //Calculare padding
    int padding;
    if (header.Width % 4 != 0)
        padding = 4 - (3 * header.Width) % 4;
    else
        padding = 0;

    int i, j;
    //Scriere pixeli in fisier
    for (i = 0; i < header.Height; i++)
    {
        for (j = 0; j < header.Width; j++)
        {
            fwrite(&matrice_pixeli[i][j], sizeof(pixelBGR), 1, f);
        }
        //Se pune paddingul
        fwrite(&zero, padding, 1, f);
    }

    fclose(f);
}

//Functie pentru modulul de recunoastere a cifrelor scrise de mana
void recunoastere(char*source,char*source_mod,char*sablon0,char*sablon1,char*sablon2,char*sablon3,char*sablon4,char*sablon5,char*sablon6,char*sablon7,char*sablon8,char*sablon9)
{
    //Declarari
    pixelBGR **matrice_source,**matrice_source_grey,**s0,**s1,**s2,**s3,**s4,**s5,**s6,**s7,**s8,**s9;
    HeaderInfo header_source,header_sablon;

    //Se citesc fisierele
    matrice_source=citire_mat(source,&header_source);
    matrice_source_grey=citire_mat(source,&header_source);
    s0=citire_mat(sablon0,&header_sablon);
    s1=citire_mat(sablon1,&header_sablon);
    s2=citire_mat(sablon2,&header_sablon);
    s3=citire_mat(sablon3,&header_sablon);
    s4=citire_mat(sablon4,&header_sablon);
    s5=citire_mat(sablon5,&header_sablon);
    s6=citire_mat(sablon6,&header_sablon);
    s7=citire_mat(sablon7,&header_sablon);
    s8=citire_mat(sablon8,&header_sablon);
    s9=citire_mat(sablon9,&header_sablon);

    //Se transforma in poze grayscale
    greyscale(matrice_source_grey,header_source);
    greyscale(s0,header_sablon);
    greyscale(s1,header_sablon);
    greyscale(s2,header_sablon);
    greyscale(s3,header_sablon);
    greyscale(s4,header_sablon);
    greyscale(s5,header_sablon);
    greyscale(s6,header_sablon);
    greyscale(s7,header_sablon);
    greyscale(s8,header_sablon);
    greyscale(s9,header_sablon);

    //Declarare si alocare vector de detectii
    detectie *D=(detectie*)malloc(100000*sizeof(detectie));

    uint32_t k=0;
    //Se genereaza detectiile pentru fiecare sablon

    printf("Se gasesc cifrele 0\n");
    matching(matrice_source_grey,s0,matrice_source,header_source,header_sablon,D,&k,0.5,CULOARE_0);
    printf("Se gasesc cifrele 1\n");
    matching(matrice_source_grey,s1,matrice_source,header_source,header_sablon,D,&k,0.499,CULOARE_1);
    printf("Se gasesc cifrele 2\n");
    matching(matrice_source_grey,s2,matrice_source,header_source,header_sablon,D,&k,0.5,CULOARE_2);
    printf("Se gasesc cifrele 3\n");
    matching(matrice_source_grey,s3,matrice_source,header_source,header_sablon,D,&k,0.5,CULOARE_3);
    printf("Se gasesc cifrele 4\n");
    matching(matrice_source_grey,s4,matrice_source,header_source,header_sablon,D,&k,0.5,CULOARE_4);
    printf("Se gasesc cifrele 5\n");
    matching(matrice_source_grey,s5,matrice_source,header_source,header_sablon,D,&k,0.5,CULOARE_5);
    printf("Se gasesc cifrele 6\n");
    matching(matrice_source_grey,s6,matrice_source,header_source,header_sablon,D,&k,0.5,CULOARE_6);
    printf("Se gasesc cifrele 7\n");
    matching(matrice_source_grey,s7,matrice_source,header_source,header_sablon,D,&k,0.5,CULOARE_7);
    printf("Se gasesc cifrele 8\n");
    matching(matrice_source_grey,s8,matrice_source,header_source,header_sablon,D,&k,0.5,CULOARE_8);
    printf("Se gasesc cifrele 9\n");
    matching(matrice_source_grey,s9,matrice_source,header_source,header_sablon,D,&k,0.5,CULOARE_9);

    free(s0);
    free(s1);
    free(s2);
    free(s3);
    free(s4);
    free(s5);
    free(s6);
    free(s7);
    free(s8);
    free(s9);

    //Se sorteaza detectiile
    qsort(D,k,sizeof(detectie),cmp);

    //Se elimina detectiile inutile
    eliminare(D,&k,header_sablon);

    //Se deseneaza conturul pentru fiecare detectie ramasa
    uint32_t i;
    for(i=0;i<k;i++)
        incadrare(matrice_source, header_source ,header_sablon.Height,header_sablon.Width,D[i].x,D[i].y,D[i].red,D[i].green,D[i].blue);

    printf("\nS-au detectat toate sabloanele\n");

    //Afisare rezultat
    afisare_mat(source_mod,header_source,matrice_source);

    free(matrice_source);
    free(matrice_source_grey);

    free(D);
}

int main()
{
    //Criptare / Decriptare
    char *numepoza_orig, *numepoza_encript,*numepoza_decript, *secretKey;
    FILE *f;
    f=fopen("fisiere.txt","r");
    numepoza_orig = (char*)malloc(20 * sizeof(char));
    numepoza_encript = (char*)malloc(20 * sizeof(char));
    numepoza_decript = (char*)malloc(20 * sizeof(char));
    secretKey = (char*)malloc(20 * sizeof(char));

    fscanf(f,"%20s", numepoza_orig);

    fscanf(f,"%20s", numepoza_encript);

    fscanf(f,"%20s", numepoza_decript);

    fscanf(f,"%20s", secretKey);

    criptare(numepoza_orig, numepoza_encript, secretKey);
    decriptare(numepoza_decript, numepoza_encript, secretKey);

    //Recunoastere pattern
    char *source=(char*)malloc(20 * sizeof(char));
    char *sablon0=(char*)malloc(20 * sizeof(char));
    char *sablon1=(char*)malloc(20 * sizeof(char));
    char *sablon2=(char*)malloc(20 * sizeof(char));
    char *sablon3=(char*)malloc(20 * sizeof(char));
    char *sablon4=(char*)malloc(20 * sizeof(char));
    char *sablon5=(char*)malloc(20 * sizeof(char));
    char *sablon6=(char*)malloc(20 * sizeof(char));
    char *sablon7=(char*)malloc(20 * sizeof(char));
    char *sablon8=(char*)malloc(20 * sizeof(char));
    char *sablon9=(char*)malloc(20 * sizeof(char));
    char *source_mod=(char*)malloc(20 * sizeof(char));

    strcpy(source,numepoza_decript);

    fscanf(f,"%20s", sablon0);
    fscanf(f,"%20s", sablon1);
    fscanf(f,"%20s", sablon2);
    fscanf(f,"%20s", sablon3);
    fscanf(f,"%20s", sablon4);
    fscanf(f,"%20s", sablon5);
    fscanf(f,"%20s", sablon6);
    fscanf(f,"%20s", sablon7);
    fscanf(f,"%20s", sablon8);
    fscanf(f,"%20s", sablon9);
    fscanf(f,"%20s", source_mod);

    recunoastere(source,source_mod,sablon0,sablon1,sablon2,sablon3,sablon4,sablon5,sablon6,sablon7,sablon8,sablon9);

    //eliberare memorie
    free(numepoza_orig);
    free(numepoza_encript);
    free(numepoza_decript);
    free(secretKey);
    free(source);
    free(sablon0);
    free(sablon1);
    free(sablon3);
    free(sablon5);
    free(sablon7);
    free(sablon9);
    free(sablon2);
    free(sablon4);
    free(sablon6);
    free(sablon8);
    return 0;
}
