#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include<math.h>

double rounded(double num, int dec){
	double pot=1;
	int temp,i;
	double arredondado;
	
	for(i=0; i<dec; i++){
		pot=pot*10;
	}
	temp=num*pot;
	arredondado=temp/pot;
	
	return arredondado;
}

// Retorna o valor em modulo
double absoluto(double x){
	if (x >= 0.0) return x;
   	else return -x;
}

// Exibe a matriz
void exibeMat(int tam, double mat[tam][tam+1]){
	int i,j;
	for(i=0;i<tam;i++){
		for(j=0;j<=tam;j++){
			if(j==0){
				printf("[ %.2lf  ",mat[i][j]);
				
			}else if(j==tam){
				printf(":  %.2lf]",mat[i][j]);
			}else{
				printf("%.2lf  ",mat[i][j]);
			}
		}
		printf("\n");
	}
}

// Limpa a matriz
void limpaMat(int tam, double mat[tam][tam+1]){
	int i,j;
	for(i=0;i<tam;i++){
		for(j=0;j<=tam;j++){
			mat[i][j] = 0;
		}
	}
}

void pivoTotal(int tam, double matriz[tam][tam+1]){
	int i,j,k,l,c;
	double mult,max,aux;
	for (k=0;k<tam-1; k++){
		max = absoluto(matriz[k][k]);
		l = c = k;
		for(i=k+1;i<tam;i++){
			for(j=k+1;j<tam;j++){
				if(absoluto(matriz[i][j]) > max){
					max = absoluto(matriz[i][j]);
					l = i;
					c = j;
				}
			}
		}
		if(l != k){
			for(j=k;j<=tam;j++){
				aux = matriz[k][j];
				matriz[k][j] = matriz[l][j];
				matriz[l][j] = aux;
			}
		}
		if(c != k){
			for(j=0;j<tam;j++){
				aux = matriz[j][k];
				matriz[j][k] = matriz[j][c];
				matriz[j][c] = aux;
			}
		}
		
		for(i=k+1;i<tam;i++){
			mult = matriz[i][k]/matriz[k][k];
			//matriz[i][k] = 0;
			for(j=k;j<=tam;j++){
				matriz[i][j] -= mult * matriz[k][j];
			}
		}	
	}
	
}

void pivotamento(int tam, double mat[tam][tam+1], int lin,int col){
	double maior=absoluto(mat[lin][col]);
	int i=0,posit=lin;
	// Encontra a posicao do maior pivo da coluna x
	for(i=lin;i<tam;i++){
		if(absoluto(mat[i][col])>maior){
			maior = absoluto(mat[i][col]);
			posit = i;
		}
	}
	// Se encontrar um maior, troca os valores
	if(posit != lin){
		double v[tam+1];
		for(i=0;i<=tam;i++){
			v[i] = mat[posit][i];				// Passa para um vetor aux
			mat[posit][i] = mat[lin][i];		// Troca
			mat[lin][i] = v[i];
		}
	}
}

bool verificaSist(int tam,double mat[tam][tam+1],int Li){
	int j=0;
	while(rounded(mat[Li][j],6) == 0 && j<tam){	// Na matriz de ordem 12, o maximo e precisao 15
		mat[Li][j] = 0;
		j++;
	}
	if(j==tam){
		if(mat[Li][tam] == 0){
			printf("\nSistema Possivel Indeterminado (SPI) \n");
		}else{
			printf("\nSistema Impossivel (SI) \n");
		}
		printf("Devido a L%d = [ ",Li);
		int i;
		for(i=0;i<=tam;i++){
			if(i==tam){
				printf(": %lf ",mat[Li][i]);
			}else{
				printf("%lf ",mat[Li][i]);
			}
		}
		printf("]\n");
		return false;
	}
	return true;
}

void GaussJordan(int tam, double L[tam][tam+1]){
	int i=0,j=0,k=0;
	bool continua=true;
	double v[tam+1],ajj=0,aij=0;
	printf("Solucao utilizando o Metodo de Gauss-Jordan e pivotamento parcial: \n");
	for(j=0;j<tam;j++){
		pivotamento(tam,L,j,j);
		ajj = L[j][j];
		for(k=0;k<=tam;k++){
			L[j][k] /= ajj;					// Lj <- Lj / ajj
			if(L[j][k] == -0){
				L[j][k] = 0;
			}
			v[k] = L[j][k];					// V <- Lj
		}
		for(i=0;i<tam && continua==true;i++){
			if(i!=j){
				aij = L[i][j];
				for(k=0;k<=tam;k++){
					L[j][k] *= aij;			// Lj <- Lj * aij
					L[i][k] -= L[j][k];		// Li <- Li * ajj
					L[j][k] = v[k];			// Lj <- V
				}
				continua = verificaSist(tam,L,i);
			}
		}
	}
	
	if(continua==true){
		printf("[A':b']:\n");
		exibeMat(tam,L);
		printf("\nSistema Possivel Determinado (SPD)\nS={( ");
		for(i=0;i<tam;i++){
			if(i==tam-1){
				printf("%.10lf",L[i][tam]);
			}else{
				printf("%.10lf ; ",L[i][tam]);
			}
		}
		printf(" )}\n\n");
	}
}

void Sassenfeld(int tam, double mat[tam][tam+1]){
	int i,j,k;
	double alp[tam],aux[tam],maior=0;
	// Limpa vetor
	for(i=0;i<tam;i++){
		alp[i] = 1;		// Vetor do valor do alpha da linha i
		aux[i] = 0;		// Vetor da soma dos coeficientes da linha i
	}
	// Calcula o beta da linha 1
	if(absoluto(mat[0][0]) != 0 && absoluto(alp[0]) != 0){			// Se o valor da diagonal principal da linha 1 for 0, atribui 0 direto
		for(i=1;i<tam;i++){
			aux[0] += absoluto(mat[0][i]);			// Senao, calcula a soma total 
		}
		alp[0]  = aux[0] / absoluto(mat[0][0]);
	}else{
		alp[0] = 0;
	}
	// Calcula das demais linhas
	for(i=1;i<tam;i++){
		if(mat[i][i]!=0){		// Se o valor da diagonal principal da linha i for diferente de 0
			for(j=0;j<tam;j++){			// For para percorrer as colunas e o vetor do alpha
				if(i!=j){
					if(absoluto(mat[i][j]) != 0 && alp[j] != 0){
						aux[i] += round(absoluto(mat[i][j]) * alp[j]);
						//aux[i] = rounded(aux[i],10);
					}
				}
			}
			if(aux[i] != 0 && absoluto(mat[i][i]) != 0){
				alp[i] = round(aux[i] / absoluto(mat[i][i]));	// Divide pelo coeficiente da diagonal principal
			}else{
				alp[i] = 0;
			}
			
		}else{
			alp[i] = 0;
		}
	}
	maior=alp[0];
	for(i=1;i<tam;i++){
		if(alp[i]>maior){
			maior=alp[i];
		}
	}
	if(maior>=1){
		printf("\nNao ha certeza que o Metodo de Gauss-Seidel convergira para a solucao do sistema, pois: ");
		printf("%.2lf > 1\n",maior);
	}else{
		printf("\no Metodo de Gauss-Seidel convergira para a solucao do sistema,pois: ");
		printf("%.2lf < 1\n",maior);
	}
}

void GaussSeidel(int tam, double mat[tam][tam+1],int k, double epsilon){
	double maiorX,Xi[tam],Xi_k[tam],soma=0,maior=0;
	int i=0,j=0,aux_k=0;
	printf("Solucao utilizando Gauss-Seidel e pivotamento completo:\n");
	pivoTotal(tam,mat);
	printf("Matriz com pivotamento completo [A':b']\n");
	exibeMat(tam,mat);
	Sassenfeld(tam,mat);
	
	printf("\nSistema x = Fx + d: \n");	
	for(i=0;i<tam;i++){			// Exibe o sistema Fx + d
		Xi[i] = Xi_k[i] = 0;
		printf("x%d= ( %.1lf ",i+1,mat[i][tam]);
		for(j=0;j<tam;j++){
			if(i!=j){
				if(mat[i][j] != 0){
					mat[i][j] *= -1;		// Troca o sinal dos coeficientes
					if(mat[i][j] > 0){		// Exibe apenas os coeficientes diferente de 0  (0x2)
						if(mat[i][j]==1){
							printf("+ x%d ",j+1);
						}else{
							printf("+ %.1lfx%d ",mat[i][j],j+1);
						}
					}else if(mat[i][j] < 0){
						if(mat[i][j] == -1){		// Nao exibe 1x2
							printf("- x%d ",j+1);			// (- x2)
						}else{
							printf("- %.1lfx%d ",fabs(mat[i][j]),j+1);	// (+ x2.1)
						}
					}		
				}				
			}
		}
		printf(") / %.2lf\n",mat[i][i]);
	}
	system("pause");
	i=j=0; aux_k = 0, maiorX=epsilon;
	for(aux_k=0;aux_k < k && maiorX >= epsilon;aux_k++){
		printf("\nK = %d\n",aux_k+1);
		for(i=0;i<tam;i++){
			soma = mat[i][tam];
			for(j=0;j<tam;j++){
				if(i!=j){
					if(mat[i][j]!=0 && Xi[j]!=0){
						soma += mat[i][j] * Xi[j];
					}
				}
			}
			if(mat[i][i] != 0 && soma != 0){
				Xi[i] = soma / mat[i][i];
			}else{
				Xi[i] = 0;
			}
			
			printf("x%d= %lf\n",i+1,Xi[i]);
		}
		
		maior = absoluto(absoluto(Xi[0]) - absoluto(Xi_k[0]));
		for(i=1;i<tam;i++){
			if(absoluto(absoluto(Xi[i]) - absoluto(Xi_k[i])) > maior){
				maior = absoluto(absoluto(Xi[i]) - absoluto(Xi_k[i]));
			}
		}
		for(j=0;j<tam;j++){
			Xi_k[j] = Xi[j];	// Passa os valores de Xi atual para que na prox iteracao sejam o Xi anterior
		}
		if(maior >= epsilon){
			printf("Maior Xk - Xk-1 = %.19lf >= %.19lf\n",maior,epsilon);
		}else{
			printf("Maior Xk - Xk-1 = %.19lf < %.19lf\n",maior,epsilon);
		}		
		maiorX = maior;
	}

	if(maiorX < epsilon){
		printf("\nCondicao de parada epsilon atingido com K=%d\n",aux_k);
		printf("%.19lf < %.19lf\n",maiorX,epsilon);
		printf("\nValores de Xi obtidos:\n");
		for(i=0;i<tam;i++){
			printf("X%d = %.19lf \n",i+1,Xi[i]);
		}
	}else{
		printf("\nO valor maximo de K foi atingido e a condicao de parada epsilon nao foi atingida \n");
		printf("%.19lf > %.19lf\n",maiorX,epsilon);
		printf("\nValores de Xi arredondados obtidos:\n");
		for(i=0;i<tam;i++){
			printf("X%d = %.19lf \n",i+1,Xi[i]);
		}
	}
}

int main(){
	FILE *arq;
	arq = fopen("Input14.txt","rt");
	if(arq == NULL){
		printf("Erro ao ler o arquivo\n");
		return 0;
	}
	char linha[250],output[250];
	char *token;
	int tam=0,i=0,j=0,k=0,valK=0;
	// Leitura da primeira linha, que contém o tamanho do sistema
	fgets(linha,250,arq);	
	tam=atoi(linha);
	double sistLin[tam][tam+1], aux[tam][tam+1],B[tam];	// Criação da matriz e do vetor com os valores de 'B'
	double epsilon=0;
	// limpar as variaveis
	limpaMat(tam,sistLin);
	limpaMat(tam,aux);
	
	while(!feof(arq)){		// Enquanto não encerrar o arquivo
		fgets(linha,250,arq);
		if(linha[0] != '\n' && linha[0] != 'O'){	// Se nao for a linha com o output ou linha vazia, entao eh a matriz
			// Numeros
			if(i<tam){		// Enquanto nao terminar de ler a matriz inteira
				token = strtok(linha," ");
				j=0;
				while(token != NULL){
					sistLin[i][j] = atof(token);		// Insere em uma matriz de float
					token = strtok(NULL," ");
					j++;
				}
				i++;
			}else{	// Senao, eh a linha com os valores independentes
				token = strtok(linha," ");
				while(token!=NULL){
					B[k] = atof(token);			// Armazena em um vetor separado
					token = strtok(NULL," ");
					k++;
				}
				// Leitura de K e Epsilon
				fgets(linha,250,arq);
				valK = atoi(linha);
				fgets(linha,250,arq);
				epsilon = atof(linha);
			}	
		}else if(linha[0] == 'O'){
			strcpy(output,linha);
		}
		// Output e linhas vazias	
	}	
	fclose(arq);
	
	// Junta os valores independentes com a matriz
	for(i=0;i<tam;i++){
		sistLin[i][tam] = B[i];
	}
	
	// Passa os valores para uma matriz auxiliar
	for(i=0;i<tam;i++){
		for(j=0;j<=tam;j++){
			aux[i][j] = sistLin[i][j];
		}
	}
	// K e Epsilon
	valK=10; epsilon=0.1;
	
	// Exibe matriz aumentada original
	printf("Matriz original [A:b]: \n");
	exibeMat(tam,sistLin);
	printf("K=%d\nEpsilon=%lf\n%s\n",valK,epsilon,output);
	printf("\n");
	
	GaussJordan(tam, sistLin);
	system("pause");
	
	GaussSeidel(tam,aux,valK,epsilon);
	
	return 0;
}
