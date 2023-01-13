/*
================================================================================
  dados.txt -s 0.3 -m 0.002 
 Paulo Alceu d´ Almeida Rezende
Projeto: Descoberta de Conhecimento em Bases de Dados em Regras de Associação Negativas
================================================================================
Implementacao do Algoritmo APRIORI NEGATIVO
================================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*------------------------------------------------------------------------------
  DEFINICAO DAS ESTRUTURAS DE DADOS
  ----------------------------------------------------------------------------*/

// tipo utilizado para a estrutura da arvore de prefixos
typedef struct t_node
{
   unsigned int sizelevel;
   struct t_item *itemlist;
} t_c;

// tipo utilizado para cada item na arvore de prefixos
typedef struct t_item
{
   unsigned int itemid,
                support,
                maxlengh;
   struct t_node *nextlevel;
   struct t_item  *nextitem;
} t_i;

// tipo utilizado para manter as transacoes na memoria
typedef struct t_trans
{
   unsigned int transid;
   unsigned transsize;
   unsigned int *itemvet;
   struct t_trans *nexttrans;
} t_t;

/*------------------------------------------------------------------------------
  DECLARACAO DE VARIAVEIS GLOBAIS
  ----------------------------------------------------------------------------*/
char vetor[1000];
char vetoraux[1000][1000];
float vetorauxsup[1000][1000];
int vetorauxconj[1000];
char filename[25];       /* string para nome de arquivo (entrada ou saida) */
int freqxy;

unsigned int numtrans, conta;      // numero de transacoes da base de dados

float percminsup,          // suporte minimo
      minMI;              // MI minimo
FILE *debfile;           // arquivo de saida

/*------------------------------------------------------------------------------
  Procedimento: Libera a memoria utilizada pela transacao
  ----------------------------------------------------------------------------*/
void freetrans(struct t_trans *auxtrans)
{
    free(auxtrans->itemvet);
    free(auxtrans);
}

/*------------------------------------------------------------------------------
  Funcao: Conta o número de transações de cada item para montar a arvore
  ----------------------------------------------------------------------------*/
unsigned int countsupport(unsigned int *vetauxitem, struct t_trans *transroot,
                        int k)
{
   struct t_trans *transaux, *transaux2;
   unsigned int support = 0, *auxitemlist, i, j;

   transaux = transroot;
   if (transaux == NULL)
      transaux2 = transaux;
   else
      transaux2 = transaux->nexttrans;
   while (transaux2 != NULL)
   {
      if (transaux2->transsize >= k)
      {
         auxitemlist = transaux2->itemvet;
         j = 0;
         for (i = 0; i < k; i++)
            while ((vetauxitem[i] != auxitemlist[j])
                    && (j < transaux2->transsize))
               j++;
         if (vetauxitem[i-1] == auxitemlist[j])
            support++;
         transaux = transaux2;
         transaux2 = transaux2->nexttrans;
      }
      else  // retira a transacao da memoria
      {
         transaux->nexttrans = transaux2->nexttrans;
         freetrans(transaux2);
         transaux2 = transaux->nexttrans;
      }
   }

   return support;
}

/*------------------------------------------------------------------------------
  Funcao: Gera os itemsets frequentes de tamanho k
  ----------------------------------------------------------------------------*/
int itemsetgen(int k, struct t_node *p_root,
                unsigned int *vetauxitem, struct t_trans *transroot,
                unsigned int minsup)
{

   struct t_item *p_aux_i, *p_aux2_i, *p_new_i, *p_auxnew_i;
   struct t_node *p_aux_n;
   unsigned int support, ok;

   ok = 0;
   if (p_root != NULL)
   {
      p_aux_i = p_root->itemlist;
      while (p_aux_i != NULL)
      {
         vetauxitem[p_root->sizelevel - 1] = p_aux_i->itemid;
         /* gera os itemsets */
         if (p_root->sizelevel == (k - 1))
         {

            // obtem cada sufixo possivel e gera os candidatos
            p_aux2_i = p_aux_i->nextitem;
            p_auxnew_i = p_aux_i->nextlevel->itemlist;
            while (p_aux2_i != NULL)
            {
               vetauxitem[k-1] = p_aux2_i->itemid;
               support = countsupport(vetauxitem, transroot, k);
               if (support >= minsup) // armazena item se tiver suporte 
               {
                  p_new_i = (struct t_item *) malloc (sizeof(struct t_item));
                  p_new_i->itemid = p_aux2_i->itemid;
                  p_new_i->support = support;
                  p_new_i->maxlengh = 1;
                  p_new_i->nextitem = NULL;
                  p_aux_n = (struct t_node *) malloc (sizeof(struct t_node));
                  p_aux_n->sizelevel = k+1;
                  p_aux_n->itemlist = NULL;
                  p_new_i->nextlevel = p_aux_n;
                  if (p_auxnew_i == NULL)
                     p_aux_i->nextlevel->itemlist = p_new_i;
                  else
                     p_auxnew_i->nextitem = p_new_i;
                  p_auxnew_i = p_new_i;
                  ok = 1;
               }
               p_aux2_i = p_aux2_i->nextitem;
            }
         }
         else
         {
            // caminha pela arvore, em profundidade
            ok = itemsetgen(k, p_aux_i->nextlevel, vetauxitem,
                                  transroot, minsup) || ok;
         }
         p_aux_i = p_aux_i->nextitem;
      }
   }
   return ok;
}

/*------------------------------------------------------------------------------
  Procedimento: Armazena todos os itemsets frequentes
  ----------------------------------------------------------------------------*/
void showitemsets(struct t_node *p_root, int k, unsigned int *vetauxitem)
{
   struct t_item *p_aux_i;
   int i;

   if (p_root != NULL)
   {
      p_aux_i = p_root->itemlist;
      while (p_aux_i != NULL)
      {
         vetauxitem[p_root->sizelevel-1] = p_aux_i->itemid;
         if (p_root->sizelevel == k)
         {
            for (i = 0; i < k; i++)               
               fprintf(debfile, "%d ", vetauxitem[i]);                
               fprintf(debfile, "(%.4f)\n", (float)p_aux_i->support/(float)numtrans);                                  
/*             fprintf(debfile, " {%d}", p_root->sizelevel); // ou i
               if (p_aux_i->nextitem == 0){
                              fprintf(debfile, " {Não}");}
               else{
                              fprintf(debfile, " {Sim}");}
               conta++;    
               fprintf(debfile, " {%d}", conta);                                                                   
               fprintf(debfile, " {%d}\n", p_aux_i->itemid);        */
         }
         else   
               showitemsets(p_aux_i->nextlevel, k, vetauxitem);                  
               p_aux_i = p_aux_i->nextitem;
      }
   }
}


/*------------------------------------------------------------------------------
  Procedimento: Varrendo arvore em profundidade 
  ----------------------------------------------------------------------------*/
void varre(struct t_node *p_root, int k, int pos,unsigned int *vetauxitem)
{

  struct t_item *p_aux_i;
  int i;
  
  p_aux_i = p_root->itemlist;
  while (p_aux_i != NULL)
  { 
         vetauxitem[p_root->sizelevel] = p_aux_i->itemid;
         vetor[p_root->sizelevel] = p_aux_i->itemid;
         conta++;         
         //if (pos==conta)
         //{
                  for (i = 1; i <= (p_root->sizelevel); i++)
                  {
                  //fprintf(debfile, " %d ",vetor[i]);
                  vetoraux[conta][i] = vetor[i];   
                  vetorauxsup[conta][i] = ((float)p_aux_i->support/(float)numtrans);                                                                                                     
                  }
         //}
         //fprintf(debfile, "- %d \n",p_aux_i->support);
         varre(p_aux_i->nextlevel, k++,pos, vetauxitem);
         p_aux_i  = p_aux_i->nextitem;         
  }
}  

/*------------------------------------------------------------------------------
  Procedimento: Montras as Regras 
  ----------------------------------------------------------------------------*/
void vetorregras(struct t_node *p_root, int k, unsigned int *vetauxitem)
{

  struct t_item *p_aux_i;
  int i;
  
  p_aux_i = p_root->itemlist;
  while (p_aux_i != NULL)
  { 
         vetauxitem[p_root->sizelevel] = p_aux_i->itemid;
         vetor[p_root->sizelevel] = p_aux_i->itemid;
         for (i = 1; i <= (p_root->sizelevel); i++)
         {
                  fprintf(debfile, "%d ",vetor[i]);               
         }
         fprintf(debfile, "- %d \n",p_aux_i->itemid);
         //fprintf(debfile, "- %d \n",p_aux_i->support);
         vetorregras(p_aux_i->nextlevel, k++, vetauxitem);
         p_aux_i  = p_aux_i->nextitem;
  }
}

/*------------------------------------------------------------------------------
  Procedimento: QuickSort (Ordenação do vetor vetorauxconj), {XuY} frequentes
  ----------------------------------------------------------------------------*/
void qksort(int ilo, int ihi) {
    int pivot;		// pivot value for partitioning array
    int ulo, uhi;	// indices at ends of unpartitioned region
    int ieq;		// least index of array entry with value equal to pivot
    int tempEntry;	// temporary entry used for swapping

    if (ilo >= ihi) {
	return;
    }
    // Select a pivot value.
    pivot = vetorauxconj[(ilo + ihi)/2];
    // Initialize ends of unpartitioned region and least index of entry
    // with value equal to pivot.
    ieq = ulo = ilo;
    uhi = ihi;
    // While the unpartitioned region is not empty, try to reduce its size.
    while (ulo <= uhi) {
	if (vetorauxconj[uhi] > pivot) {
	    // Here, we can reduce the size of the unpartitioned region and
	    // try again.
	    uhi--;
	} else {
	    // Here, A[uhi] <= pivot, so swap entries at indices ulo and
	    // uhi.
	    tempEntry = vetorauxconj[ulo];
	    vetorauxconj[ulo] = vetorauxconj[uhi];
	    vetorauxconj[uhi] = tempEntry;
	    // After the swap, A[ulo] <= pivot.
	    if (vetorauxconj[ulo] < pivot) {
		// Swap entries at indices ieq and ulo.
		tempEntry = vetorauxconj[ieq];
		vetorauxconj[ieq] = vetorauxconj[ulo];
		vetorauxconj[ulo] = tempEntry;
		// After the swap, A[ieq] < pivot, so we need to change
		// ieq.
		ieq++;
		// We also need to change ulo, but we also need to do
		// that when A[ulo] = pivot, so we do it after this if
		// statement.
	    }
	    // Once again, we can reduce the size of the unpartitioned
	    // region and try again.
	    ulo++;
	}
    }
    // Now, all entries from index ilo to ieq - 1 are less than the pivot
    // and all entries from index uhi to ihi + 1 are greater than the
    // pivot.  So we have two regions of the array that can be sorted
    // recursively to put all of the entries in order.
    qksort(ilo, ieq - 1);
    qksort(uhi + 1, ihi);
}
/*------------------------------------------------------------------------------
  Funcao: Calcula o MI de uma regra Negativa
  ----------------------------------------------------------------------------*/
float MI_calc(float sup_xy, float antec_sup, float conseq_sup)
{
   float aux;
   aux =  sup_xy -(antec_sup * conseq_sup);
   if (aux == 0)
      return 0;
   else
      return (aux / antec_sup);
}

/*------------------------------------------------------------------------------
  Função: para saber se temos {XuY} frequentes
  ----------------------------------------------------------------------------*/
int verificaconjunto(int j, int l, int i)
{
int k, item;
int m;
int cont,cont2;
int ok;
int a; //para varrer o vetor para comparar
int b; //para varrer o vetor para comparar
char c;

   
ok =0;
cont = 0;
     //fprintf(debfile, "|");
     
     // Gera o vetor conjunto {XuY}
     for (k = 1; k <= (i - 2); k++)
     {
          if (vetoraux[j][k] != 0)
          {
          cont++;
          //fprintf(debfile, "%d", vetoraux[j][k]);
          vetorauxconj[cont] = vetoraux[j][k];
          }
     }
         
     for (m = 1; m <= (i - 2); m++)
     {
           if (vetoraux[l][m] != 0)
           {
           cont++;           
           vetorauxconj[cont] = vetoraux[l][m];           
           }
     }  
     
     qksort(0,cont); //ordenando o vetor conjunto {XuY} usando o quicksort
  

 /*   fprintf(debfile, "-");

      for (a = 1; a <= conta; a++)
      {
           fprintf(debfile, "|");
           for (b = 1; b <= cont; b++)
           {
               //if (vetoraux[a][b] !=0){
               fprintf(debfile, "%d", vetoraux[a][b]);//}
           }
       }*/
        
       
 /*for (b = 1; b <= cont; b++)
 {
   if (vetoraux[a][b] !=0){
               fprintf(debfile, "%d", vetoraux[a][b]);//}
 }

*/


   //Veficando a frequencia de {XuY} na base de dados para calcular o MI.
   FILE   *bd_in2;      /* arquivo de entrada da base de dados   */
   if ((bd_in2 = fopen(filename, "rb"))==NULL)
   {
      printf("\nERRO: Nao foi possivel abrir o arquivo %s.\n", filename);
      exit(1);
   }
              
   FILE   *bd_in3;      /* arquivo de entrada da base de dados   */
   if ((bd_in3 = fopen(filename, "rb"))==NULL)
   {
      printf("\nERRO: Nao foi possivel abrir o arquivo %s.\n", filename);
      exit(1);
   }  





   freqxy = 0;                  
  /* Le os itens do arquivo de transacoes e armazena-os na memoria */
   fseek(bd_in2, 0, 0);
   c = getc(bd_in2);

   while (!feof(bd_in2))
   {
      while (c != '\r')
      {
	     item = 0;
         while ( (c >= '0') && (c <= '9'))
	     {
            fscanf(bd_in3, "%d", &item);                          
            //fprintf(debfile, "!%c!", c);   
            for (m = 1; m <= cont; m++)
            {   
               if (vetorauxconj[m] ==  item)
               {
   
                    cont2++;                                   
                    if (cont2 == cont)
                    {
                      //fprintf(debfile, "ok");              
                      freqxy++;
                      break;
                    }             
             
               }
            }
            
            
            c = getc(bd_in2);
         }
	     if (c == ' ') c = getc(bd_in2);
      }
      if (c == '\r')
      {
            //fprintf(debfile, "$");     
            cont2 = 0;         
      }
      c = getc(bd_in2);
      if (!feof(bd_in2)) c = getc(bd_in2);
   }                                  
                                                                            
   fclose(bd_in2);                                                                  
                                                                          
                                                                                  
                                                                    
       
   //fprintf(debfile, "!%d!", freqxy);          
          
   //Apenas mostrando o {XuY}.
   fprintf(debfile, " - XuY={");  
     for (m = 1; m <= cont; m++)
     {
         fprintf(debfile, "%d", vetorauxconj[m]);
     }
     fprintf(debfile, "}");  
                                          
                                                                          
    // Compara o vetor do conjunto {XuY} com os vetor de itens frequentes 
    for (a = 1; a <= conta; a++)
      {
           for (b = 1; b <= cont; b++)
           {
               if (vetorauxconj[b] == vetoraux[a][b])
               {
               ok++;
               }
               else
               {
               ok--;
               }
           }
           if (ok == cont)
           {
                      return 0;
           }
           ok = 0;
       }     

     
     //Apenas mostrando o {XuY}.
 /*  fprintf(debfile, " - XuY={");  
     for (m = 1; m <= cont; m++)
     {
         fprintf(debfile, "%d", vetorauxconj[m]);
     }
     fprintf(debfile, "}");  */
     
     return 1;      
     
     
} 

void showregras(int i, float MImin)
{
   int ok,j,l,k,m;
   float antec,conc; // Para calculo de MI
   ok = 1; // variavel para saber se temos itens no antecedente que
           // existe no consequente

   // Gera todas as regras com suas condições... 
   for (j = 1; j <= conta; j++)
   {
     fprintf(debfile, "\n");   
     for (l = 1; l <= conta; l++)
     {
     fprintf(debfile, "\n");   
                for (k = 1; k <= (i - 2); k++)
                {                      
                  for (m = 1; m <= (i - 2); m++)
                  {
                      if ((vetoraux[j][k] == vetoraux[l][m]) && (vetoraux[j][k] != 0))
                      {
                      ok = 0;
                      }
                  }
                  if (vetoraux[j][k] != 0){ 
                     fprintf(debfile, " %d ", vetoraux[j][k]);
                     antec = vetorauxsup[j][k];  // armasena o sup do antecedente
                  }
                  if (k == (i - 2)) // fim da regra mostrada...
                  {
                       fprintf(debfile, " -/>");            
                       for (m = 1; m <= (i - 2); m++)
                       {
                          if (vetoraux[l][m] != 0)
                          {
                             fprintf(debfile, " %d", vetoraux[l][m]);   
                             conc = vetorauxsup[l][m];    // armasena o sup do consequente                     
                          }
                       }  
                       if (verificaconjunto(j,l,i) != 1) // função para saber se temos {XuY} frequentes
                       {
                           fprintf(debfile, " - temos {XuY} frequentes");
                       }
                       else if (MI_calc(freqxy,antec,conc)<= MImin)                      
                       {
                          fprintf(debfile, " - MI não suportado");
                       }                       
                       else if (ok!=1)                      
                       {
                          fprintf(debfile, " - {X} tem em {Y}");
                       }          
                       else
                       {
                            fprintf(debfile, " - X(%.4f) Y(%.4f) ESP(%.4f) MI(%.4f)",antec,conc,antec*conc, MI_calc(freqxy,antec,conc)); //Mostra Valores
                       }              
                       ok = 1;
                  }
              }
     }
   }
}


/*------------------------------------------------------------------------------
  Procedimento: Exibe os parametros do programa
  ----------------------------------------------------------------------------*/
void showparam()
{
   printf("\nSINTAXE: apriori_v1 <arquivo de transacoes> [opcoes]\n\n");
   printf("Opcoes: -o <arquivo de saida>  => nome do arquivo de saida\n");
   printf("        -s <minsup>            => valor do suporte minimo\n");
   printf("        -m <minMI>           => valor do MI minimo\n\n");
   printf("Exemplo: apriori_v1 arqtrans.txt regras.txt -s 0.04 -m 0.1\n\n");
   printf("Os valores default utilizados para as opcoes nao informadas sao:\n");
   printf("        arquivo de saida = [arquivo de transacoes].out\n");
   printf("        minsup =  0.03\n");
   printf("        minMI = 0.0002\n\n");
   printf("Caso seja fornecido o MI, as regras obtidas irao\n");
   printf("atender a medida.\n");
}

/*------------------------------------------------------------------------------
  PROGRAMA PRINCIPAL -----------------------------------------------------------
  ----------------------------------------------------------------------------*/

void main (int argc, char **argv) {

   /* Declaracao de variaveis */
   FILE                 *bd_in;      /* arquivo de entrada da base de dados   */
   struct t_node *p_root, *p_aux_n;
   struct t_item *p_new_i, *p_aux_i;
   struct t_trans *transroot, *transaux1, *transaux2;

   char fileout[25],        /* string para nome de arquivo de saida           */
        argval[4];
   unsigned int
        numitems,           /* numero de itens                                */
        maxitems,
        i, j, l, k, m, ok,  /* variavel auxiliar para contador                */
        totallarge,         /* total de itemsets frequentes                   */
        minsup;             /* suporte minimo - em no. de transacoes          */
   unsigned int item, auxitem;
   char c;
   unsigned int *vetconta, *vetordena, *vetauxitem, *vetitemlist;

   // Responsável para mostrar a Hora e a Data atual.
   char dateStr [9]; //Data e Hora
   char timeStr [9];
   _strdate( dateStr);
   _strtime( timeStr );

   /* --------------------------------------------------------------------------
      Trata os parametros de entrada do algoritmo
      ------------------------------------------------------------------------*/
   /* Inicializa variaveis */
   numitems = 0;
   numtrans = 0;
   maxitems = 0;
   /* Verifica os parametros fornecidos */
   if ((argc < 2) || (argc > 8))
   {
      printf("\nERRO: Numero incorreto de parametros.\n");
      showparam();
      exit(1);
   }
   else
   {
      /* Obtem nome do arquivo de transacoes */
      sprintf(filename, "%s", argv[1]);
      /* Inicializa parametros */
      percminsup = 0.03;
      minMI = 0.0002;
      
      sprintf(fileout, "%s.out", filename);
      /* obtem os parametros */
      i = 2;
      while (i < argc) {
         strcpy(argval, argv[i++]);
         switch (argval[1]) {
            case 'o': // arquivo de saida
               strcpy(fileout, argv[i++]);
               break;
            case 's': // valor dos suporte minimo
               percminsup = atof(argv[i++]);
               break;
            case 'm': // valor do MI minimo
               minMI = atof(argv[i++]);
               break;               
            default:
               printf("\nERRO: Parametro %d incorreto.\n", argval);
               showparam();
               exit(1);
         }
      }
   }
   /* Abre o arquivo de saida */
   debfile = fopen(fileout, "w+");

   /* --------------------------------------------------------------------------
      Le o arquivo de transacoes
      ------------------------------------------------------------------------*/
   /* Abre o arquivo com a base de transacoes */
   if ((bd_in = fopen(filename, "rb"))==NULL)
   {
      printf("\nERRO: Nao foi possivel abrir o arquivo %s.\n", filename);
      exit(1);
   }
   /* Conta o numero de itens do arquivo */
   while (!feof(bd_in))
   {
      c = getc(bd_in);
      fscanf(bd_in, "%d", &item);
      numitems++;
       if (item > maxitems) maxitems = item;
   }
   /* aloca vetor auxiliar para os items de uma trasacao */
   vetauxitem = (unsigned int *) malloc (sizeof(unsigned int) * maxitems);
   /* aloca vetor para contar itemsets de tamanho 1 */
   vetconta = (unsigned int *) malloc ((sizeof(unsigned int) * maxitems) + 1);
   for (i = 0; i <= maxitems; i++) vetconta[i] = 0;
   /* aloca a raiz da lista de transacoes */
   transroot = (struct t_trans *) malloc (sizeof(struct t_trans));
   transroot->transid = 0;
   transroot->transsize = 0;
   transroot->itemvet = NULL;
   transroot->nexttrans = NULL;
   transaux1 = transroot;
   /* Le os itens do arquivo de transacoes e armazena-os na memoria */
   fseek(bd_in, 0, 0);
   c = getc(bd_in);
   while (!feof(bd_in))
   {
      auxitem = 0;
      while (c != '\r')
      {
	     item = 0;
         while ( (c >= '0') && (c <= '9'))
	      {
            item *= 10;
            item += (int) c - (int) '0';
            c = getc(bd_in);
         }
         vetauxitem[auxitem++] = item;
	     vetconta[item]++;
	      if (c == ' ') c = getc(bd_in);
      }
      numtrans++;
      /* cria a transacao na memoria */
      if (auxitem > 1)
      {
         vetitemlist = (unsigned int *) malloc(sizeof(unsigned int) * auxitem);
         for (i = 0; i < auxitem; i++)
            vetitemlist[i] = vetauxitem[i];
         transaux2 = (struct t_trans *) malloc (sizeof(struct t_trans));
         transaux2->transid = numtrans;
         transaux2->transsize = auxitem;
         transaux2->itemvet = vetitemlist;
         transaux2->nexttrans = NULL;
         transaux1->nexttrans = transaux2;
         transaux1 = transaux2;
      }
      c = getc(bd_in);
      if (!feof(bd_in)) c = getc(bd_in);
   }
   /* Calcula o valor do suporte minimo em transacoes */
   minsup = (int) (numtrans * percminsup);

   /*--------------------------------------------------------------------------
     Exibe informacoes sobre a execucao
     -------------------------------------------------------------------------*/
   printf("\n\n");
   printf("Monografia - Paulo Alceu d´ Almeida Rezende\n");   
   printf("Agradecimentos: Programa de Iniciacao Cientifica e\r");   
   printf("\n-----------------------------------------------------\n");
   printf("Algoritmo APRIORI NEGATIVO (implementacao versao 1.0)\n");
   printf("-----------------------------------------------------\n");
   printf("** Parametros **\n");
   printf("\tArquivo de transacoes: %s\n", filename);
   printf("\tNumero de transacoes : %d\n", numtrans);
   printf("\tNumero de itens: %d\n", numitems);
   printf("\tMedia de itens por transacao: %.2f\n", (float)numitems/(float)numtrans);
   printf("\tMaior identificador de item: %d\n", maxitems);
   printf("\tSuporte minimo : %.4f (%d transacoes)\n", percminsup, minsup);
   printf("\tMI minimo : %.4f\n\n", minMI);   
   printf("** Execucao iniciada **\n\n");

   /*--------------------------------------------------------------------------
     Grava cabecalho do arquivo de saida
     -------------------------------------------------------------------------*/
   fprintf(debfile, "Paulo Alceu d´ Almeida Rezende\n");      
   fprintf(debfile, "\n-----------------------------------------------------\n");
   fprintf(debfile, "Algoritmo APRIORI NEGATIVO (implementacao versao 1.0)\n");
   fprintf(debfile, "-----------------------------------------------------\n");
   fprintf(debfile, "Executado em: Data: %s Hora: %s\n", dateStr,timeStr);   
   fprintf(debfile, "** Parametros **\n");
   fprintf(debfile, "   Arquivo de transacoes: %s\n", filename);
   fprintf(debfile, "   Numero de transacoes : %d\n", numtrans);
   fprintf(debfile, "   Numero de itens: %d\n", numitems);
   fprintf(debfile, "   Media de itens por transacao: %.2f\n", (float)numitems/(float)numtrans);
   fprintf(debfile, "   Maior identificador de item: %d\n", maxitems);
   fprintf(debfile, "   Suporte minimo : %.4f (%d transacoes)\n", percminsup, minsup);
   fprintf(debfile, "   MI minimo : %.4f\n\n", minMI);   
   fprintf(debfile, "** Resultados apos Execucao **\n\n");

   /*--------------------------------------------------------------------------
     Obtem os itemsets frequentes de tamanho 1
     -------------------------------------------------------------------------*/
   _strtime( timeStr );     
   printf("\t> Gerando os itemsets: Hora: %s\n",timeStr);
   printf("\t     .. de tamanho 1\n");
   totallarge = 0;
   if (minsup > 0)
      for (i = 0; i <= maxitems; i++)
         if (vetconta[i] >= minsup) 
	 { 
            totallarge++;
            //fprintf(debfile, "%d (%d)\n", i, vetconta[i]);
	 }
   /* cria vetor auxiliar para ordenacao */
   vetordena = (unsigned int *) malloc (sizeof(unsigned int) * totallarge);
   /* alimenta vetor auxiliar ordenado */
   j = 0;
   if (minsup > 0)
      for (i = 0; i <= maxitems; i++)
         if (vetconta[i] >= minsup)
            vetordena[j++] = i;  
   // cria a raiz da arvore de prefixos
   p_root = (struct t_node *) malloc (sizeof(struct t_node));
   p_root->sizelevel = 1;
   p_root->itemlist = NULL;
   // gera o primeiro nivel da arvore => os 1-itemsets
   for (i = 0; i < totallarge; i++)
   {
      p_new_i = (struct t_item *) malloc (sizeof(struct t_item));
      p_new_i->itemid = vetordena[i];
      p_new_i->support = vetconta[vetordena[i]];
      p_new_i->maxlengh = 1;
      p_new_i->nextitem = NULL;
      p_aux_n = (struct t_node *) malloc (sizeof(struct t_node));
      p_aux_n->sizelevel = 2;
      p_aux_n->itemlist = NULL;
      p_new_i->nextlevel = p_aux_n;
      if (i == 0)
         p_root->itemlist = p_new_i;
      else
         p_aux_i->nextitem = p_new_i;
      p_aux_i = p_new_i;
   }

   /*--------------------------------------------------------------------------
     Obtem os itemsets frequentes de tamanho maior do que 1
     -------------------------------------------------------------------------*/
   i = 2;
   do
   {
      ok = 0;
      printf("\t     .. de tamanho %d\n", i);
      ok = ok || itemsetgen(i++, p_root, vetauxitem, transroot, minsup);
   }
   while (ok);

   /*--------------------------------------------------------------------------
     Grava os itemsets frequentes no arquivo de saida
     -------------------------------------------------------------------------*/
   // i-2 corresponde ao tamanho maximo dos itemsets frequentes
   fprintf(debfile, "--> Itemsets frequentes (Suporte) <--\n");
   //fprintf(debfile, "--> Itemsets frequentes (Suporte) (Tamanho) {Usa no prox conj} {nº} {Ultimo item}<--\n");   
   for (j = 1; j <= (i - 2); j++)
   {
      fprintf(debfile, "   %d-itemsets   \n", j);
      showitemsets(p_root, j, vetauxitem);   
   }

   /*--------------------------------------------------------------------------
     Varrendo Arvore em Profundidade
     -------------------------------------------------------------------------*/
   _strtime( timeStr );     
   fprintf(debfile, "\n**  Varrendo Arvore em Profundidade ** Hora: %s\n", timeStr);  
   conta = 0;
   varre(p_root, j, 5, vetauxitem);
   vetorregras(p_root, j, vetauxitem);  // Preenche um vetor com os items frequentes
   
   
   /*--------------------------------------------------------------------------
     Gerando Regras
     -------------------------------------------------------------------------*/

   _strtime( timeStr );
   fprintf(debfile, "\n**  Gerando Regras ** Hora: %s\n", timeStr); 

   showregras(i,minMI);
      
   /*--------------------------------------------------------------------------
     Fecha os arquivos de trabalho e libera a memoria alocada
     -------------------------------------------------------------------------*/
   _strtime( timeStr );   
   fprintf(debfile, "\n\n** Finalizado em: ** Hora: %s\n", timeStr);
     
   fclose(debfile);
   fclose(bd_in);
   

} // FIM DO PROGRAMA
