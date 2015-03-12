#include <conio.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <locale.h>
using namespace std;
//еепааоарпблблллрп
double** A;
double* F;
int** KrUsel; // ћассив граничных узлов
double *H[7], *Xk, *Xk1;		
double eps = 1e-15;
int maxiter=2000,n,m,l;

struct Area
{
    double x1, y1, x2, y2, x3, y3, x4, y4;//координаты точек
    int nx1, nx2,nx3, ny1, ny2;           //число разбиений(отрезков)
    double hx01, hx02, hx03, hy01, hy02;//начальные шаги по x и y
};

Area area;		//область area
double** points; //массив с узлами
int pointsN;	//ќбщее количество узлов


bool input()	//ввод данных
{	
	FILE* in;
	in=fopen("input.txt","r");
	fscanf(in,"%lf %lf", &area.x1, &area.y1);
	fscanf(in,"%lf %lf", &area.x2, &area.y2);
	fscanf(in,"%lf %lf", &area.x3, &area.y3);
	fscanf(in, "%lf %lf", &area.x4, &area.y4);
	fscanf(in, "%d", &area.nx1);
	fscanf(in, "%d", &area.nx2);
	fscanf(in, "%d", &area.nx3);
	fscanf(in, "%d", &area.ny1);
	fscanf(in, "%d", &area.ny2);
	if (area.x1 < area.x3 &&	area.x2 < area.x3 && area.x3 < area.x4 &&
		area.y1 < area.y2 && area.y2 == area.y3 &&	area.y3 < area.y4 &&
		area.nx1 > 0 &&  area.nx2 > 0 && area.nx3 > 0 &&
		area.ny1 > 0 && area.ny2 > 0)			
		return true;
	else 
	{
		return false;
	}
}
//исправила, работает
//генераци€ сетки
void genSet()
{
	int i, j;
	double x, y, hx, hy;
	//вычисл€ем количество узлов
	pointsN = (area.nx2 + 1) * (area.ny1 + 1)	+ (area.nx1 + area.nx2 + area.nx3 + 1) * (area.ny2 + 1) - (area.nx2 + 1);
	points = new double*[pointsN];		
	for (i = 0; i < pointsN; i++)		//создаЄм массив с координатами узлов
		points[i] = new double[2];
	//вычисление длины шага
	area.hx01 = (area.x1 - area.x2) / area.nx1;
	area.hx02 = (area.x3 - area.x1) / area.nx2;
	area.hx03 = (area.x4 - area.x3) / area.nx3;
	area.hy01 = (area.y2 - area.y1) / area.ny1;
	area.hy02 = (area.y4 - area.y3) / area.ny2;
	////// расставл€ем координаты узлов 
	x = area.x1;
	y = area.y1;
	i = 0;
	while (y < area.y2)
	{
		while (x <= area.x3)
		{
			points[i][0] = x;
			points[i][1] = y;
			x += area.hx02;
			i++;
		}
		x = area.x1;
		y += area.hy01;
	}
	y = area.y2;
	while (y <= area.y4)
	{
		x = area.x2;
		while (x < area.x1)
		{
			points[i][0] = x;
			points[i][1] = y;
			x += area.hx01;
			i++;
		}
		x = area.x1;
		while (x < area.x3)
		{
			points[i][0] = x;
			points[i][1] = y;
			x += area.hx02;
			i++;
		}
		x = area.x3;
		while (x <= area.x4)
		{
			points[i][0] = x;
			points[i][1] = y;
			x += area.hx03;
			i++;
		}
		y += area.hy02;
	}

	F = new double[pointsN];
	A = new double*[pointsN];
	for (int i = 0; i < pointsN; i++)
		A[i] = new double[7];

	Xk = new double[pointsN];
	Xk1 = new double[pointsN];
	for (i = 0; i<pointsN; i++)
		Xk[i] = 0.;
	for (int i = 0; i<pointsN; i++)
		Xk1[i] = 0.;
	n = pointsN;
	m = area.nx1 - 1;
	l = area.nx2 - 1;
}

void output_grid()//	вывод полученной сетки на экран
{
	int i;
	cout << " оличество точек:" << pointsN << endl;
	for (i = 0; i < pointsN; i++)
		cout << "(" << points[i][0] << ", " << points[i][1] << ")" << endl;
}
/*
bool input_matrix()
{
    int i, t;
    FILE *in = fopen("grid.txt","r");
    if (!in) return false;

    fscanf(in,"%d",&pointsN);
    points = new double*[pointsN];
    F = new double[pointsN];
    A = new double*[pointsN];
 
	//создание массива
for(i=0; i<pointsN; i++)
    {
        points[i] = new double[2];
        A[i] = new double[7];
    }
	///////////////////////////////////////
    for(i=0; i<pointsN; i++)
        fscanf(in,"%lf %lf",&points[i][0],&points[i][1]);
    fclose(in);

    return true;
}*/

double U(double x, double y)
{
	//return x*x+y*y;
	//return 1.;
	return x+y;
}

double gamma(double x, double y)
{
    return 1.0;
}

double f(double x, double y)
{
	return x+y;
	//return  1.;
	//return -4.0+x*x+y*y;
}

//обнуление строки в массиве
void obnulenie(int i)
{
    for(int j=0; j<7; j++)
        A[i][j] = 0.0;
}

//исправила, вроде работает
// ѕоиск узла среди граничных ////////
int kr_poisk(int t)
{
	int i, j, k;
	for (i = 0; i<8; i++){
		switch (i){
		case 0: k = area.nx2 + 1; break;
		case 1: k = area.ny1 + 1; break;
		case 2: k = area.nx3 + 1; break;
		case 3: k = area.ny2 + 1; break;
		case 4: k = area.nx1 + area.nx2 +area.nx3 + 1; break;
		case 5: k = area.ny1 + 1; break;
		case 6: k = area.nx1 + 1; break;
		case 7: k = area.ny2 + 1; break;
		}
		for (j = 0; j<k; j++){
			if (KrUsel[i][j] == t)
				return i;
		}
	}
	return (-1);
}

//исправила, работает
//создание массива граничных узлов
void create_kr_usel()
{
	int i, t, k;
	//t - количество узлов
	//выделили пам€ть под 6 граней
	KrUsel = new int*[8];
	//выделили пам€ть под каждую грань в зависимости от кол-ва узлов в грани
	KrUsel[0] = new int[area.nx2 + 1];
	KrUsel[1] = new int[area.ny1 + 1];
	KrUsel[2] = new int[area.nx3 + 1];
	KrUsel[3] = new int[area.ny2 + 1];
	KrUsel[4] = new int[area.nx1 + area.nx2 + area.nx3 + 1];
	KrUsel[5] = new int[area.ny2 + 1];
	KrUsel[6] = new int[area.nx1 + 1];
	KrUsel[7] = new int[area.ny1 + 1];
	//перва€ грань
	t = area.nx2 + 1;
	for (i = 0; i < t; i++)
		KrUsel[0][i] = i;

	//втора€ грань
	KrUsel[1][0] = KrUsel[0][area.nx2];
	t = area.ny1;
	for (i = 1; i<t; i++)
		KrUsel[1][i] = KrUsel[1][i - 1] + (area.nx2 + 1);
	KrUsel[1][t] = KrUsel[1][t - 1] + area.nx1 + area.nx2 + 1;

	//треть€ грань
	t = area.nx3 + 1;
	KrUsel[2][0] = KrUsel[1][area.ny1];
	for (i = 1; i<t; i++)
		KrUsel[2][i] = KrUsel[2][i - 1] + 1;
	//четверта€ грань
	t = area.ny2 + 1;
	KrUsel[3][0] = KrUsel[2][area.nx3];
	for (i = 1; i<t; i++)
		KrUsel[3][i] = KrUsel[3][i - 1] + (area.nx1 + area.nx2+area.nx3 + 1);
	//п€та€ грань
	t = area.nx1 + area.nx2 +area.nx3 + 1;
	KrUsel[4][t - 1] = KrUsel[3][area.ny2];
	for (i = t - 2; i >= 0; i--)
		KrUsel[4][i] = KrUsel[4][i + 1] - 1;
	//шеста€ грань
	KrUsel[5][0] = KrUsel[4][0];
	t = area.ny2 + 1;
	for (i = 1; i<t; i++)
		KrUsel[5][i] = KrUsel[5][i - 1] - (area.nx1 + area.nx2 + area.nx3 + 1);
	//седьма€ грань
	t = area.nx1 + 1;
	KrUsel[6][0] = KrUsel[5][area.ny2];
	for (i = 1; i < t; i++)
		KrUsel[6][i] = KrUsel[6][i - 1] + 1;

	//¬осьма€ грань
	t = area.ny1;
	KrUsel[7][0] = KrUsel[6][area.ny1];
	for (i = 1; i < t; i++)
		KrUsel[7][i] = KrUsel[7][i - 1] - (area.nx1 + area.ny2 + 1);
	KrUsel[7][t] = KrUsel[0][0];
	
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 7; j++)
			cout << KrUsel[i][j] << "\t";
		cout << endl;
	}
}

//¬роде работает, но возможно немного неправильно
//обработка регул€рных узлов
void reg_usel(int i)
{
    int t;
    double hx, hx1, hy, hy1;

    //вычислили h(i) и h(i-1) по ’
    hx = points[i+1][0] - points[i][0];	//h(i)
    hx1 = points[i][0] - points[i-1][0];	//h(i-1)

    A[i][3] = -1.0*(2.0/(hx*hx1) + gamma(points[i][0],points[i][1]));//вычисл€ем U(i,j)
    A[i][2] = 2.0/(hx1*(hx+hx1));//вычисл€ем U(i-1,j)
    A[i][4] = 2.0/(hx*(hx+hx1));//вычисл€ем U(i+1,j)

    //вычислили h(i) и h(i-1) по Y
	t = (area.nx2+1)*(area.ny1+1);
		//если узел находитс€ в нижнем пр€моугольнике
		if(i<(t-area.nx2-1)){
				hy = points[i + (area.nx2+1)][1] - points[i][1];
                hy1 = points[i][1] - points[i - (area.nx2+1)][1];
                A[i][1] = 2.0/(hy1*(hy+hy1));//вычисл€ем Ui,j-1
				A[i][5] = 2.0/(hy*(hy+hy1));}//вычисл€ем Ui,j+1
		//если узел на стыке двух пр€моугольников
		if(i<(t+area.nx1) && i>=(t-area.nx2-1)){
			hy = points[i + area.nx1 + area.nx2 + area.nx3 + 1][1] - points[i][1];
			hy1 = points[i][1] - points[i - (area.nx1 + area.nx2 + 1)][1];
            A[i][1] = 2.0/(hy1*(hy+hy1));//вычисл€ем Ui,j-1
            A[i][6] = 2.0/(hy*(hy+hy1));//вычисл€ем Ui,j+1
        }
		//верхний пр€моугольник
		if(i>t){
			hy = points[i+area.nx1+area.nx2+area.nx3+1][1] - points[i][1];
			hy1 = points[i][1] - points[i-(area.nx1+area.nx2+area.nx3+1)][1];
            A[i][0] = 2.0/(hy1*(hy+hy1));//вычисл€ем Ui,j-1
            A[i][6] = 2.0/(hy*(hy+hy1));//вычисл€ем Ui,j+1
            } 

    A[i][3] += -2.0/(hy*hy1);//вычисл€ем Ui,j y-€ добавка
	F[i] = -f(points[i][0],points[i][1]);
}
//первые краевые
void FirstBorderType(int i)
{
    A[i][3] = 1.0;
    F[i] = U(points[i][0],points[i][1]);
}


//создание —Ћј”
void create_SLAU()
{
    int i,p;
    create_kr_usel();   //создаем массив с вершинами, которые на границе
    for(i=0; i<pointsN; i++)
	{ 
        obnulenie(i); 
        F[i] = 0.;
        if((p=kr_poisk(i))==-1) //если i-€ вершина не на границе
            reg_usel(i);      //то обрабатываем ее как регул€рную
        else
            FirstBorderType(i); //иначе - как краевую
	}
	for (int i = 0; i < pointsN; i++)
	{
		for (int j = 0; j < 7; j++)
			cout << A[i][j] << "\t";
		cout << endl;
	}
	for (int i = 0; i < pointsN; i++)
		cout << F[i]<<endl;
}
//вывод матрицы
void output_matrix()
{
    int i,j;
	
    FILE *out = fopen("in.txt","w");
	
	fprintf(out,"%d %d %d \n ", pointsN, area.nx1-1, area.nx2-1);

	for (j=0;j<7;j++)
	{for(i=0; i<pointsN; i++)
		fprintf(out,"%lf  ",A[i][j]);
	fprintf(out,"\n");}
	for(i=0; i<pointsN; i++)
		fprintf(out,"%lf\n  ",F[i]);
	fclose(out);
}
void New_input_matrix_1()
{
	int i,j;
	FILE *f = fopen("in.txt","r");
	fscanf(f,"%d%d%d",&n,&m,&l);
	
	Xk=new double[n];
	Xk1=new double[n];
	
	for(i=0; i<7; i++)
	{
		H[i]=new double[n];
		for(j=0; j<n; j++)
			{fscanf(f,"%lf",&H[i][j]);
			 printf("%lf\n ",H[i][j]);}
		printf("\n\n ");
	}
	for(int i=0; i<n; i++)
		fscanf(f,"%lf",&F[i]);
	for(i=0; i<n; i++)
		Xk[i]=0.;
	
	for(int i=0; i<n; i++)
		Xk1[i]=0.;

	fclose(f);
}

void iter(double w, double *Xk, double *Xk1)
{
	 double s=0.;
	for(int i=0; i<n; i++)
		{
		 s=0.;
		s=H[3][i]*Xk[i];
		if(i>0)
			s+=H[2][i]*Xk[i-1];
		if(i<n-1)
			s+=H[4][i]*Xk1[i+1];
		if(i+2+m<n)
			s+=H[5][i]*Xk1[i+2+m];
		if(i-2-m>=0)
			s+=H[1][i]*Xk[i-2-m];
		if(i+3+l+m<n)
			s+=H[6][i]*Xk[i+3+l+m];
		if(i-3-l-m>=0)
			s+=H[0][i]*Xk[i-3-l-m];
		Xk1[i]=Xk[i]+(w/H[3][i])*(F[i]-s);
		}
}
		
double nevjazka()
{
	double res=0., res1=0.;
	for(int i=0; i<n; i++)
	{
		double s=F[i];
		s-=H[3][i]*Xk[i];
		if(i>0)
			s-=H[2][i]*Xk[i-1];
		if(i<n-1)
			s-=H[4][i]*Xk[i+1];
		if(i+2+m<n)
			s-=H[5][i]*Xk[i+2+m];
		if(i-2-m>=0)
			s-=H[1][i]*Xk[i-2-m];
		if(i+3+l+m<n)
			s-=H[6][i]*Xk[i+3+l+m];
		if(i-3-l-m>=0)
			s-=H[0][i]*Xk[i-3-l-m];
		res+=s*s;
		res1+=F[i]*F[i];
	}

	res=sqrt(res);
	res1=sqrt(res1);
	res=res/res1;
	return res;
}


void Jakoby(double w)
{
	for(int it=1; it<=maxiter; it++)
	{
		iter(w, Xk, Xk1);
		for(int i=0; i<n; i++)
		Xk[i]=Xk1[i];
		double nev=nevjazka();
		printf("%d: %lf\n",it,nev);
		if(nev<=eps) 
			break;
	}

}

void WriteRes()
{double norm=0.0, norm1=0.0;
	FILE *fp = fopen("res.txt", "w");
	for(int i=0; i<n; i++)
		fprintf(fp, "%le\n", Xk[i]);
	fprintf(fp, "\n\n");
	for(int i=0; i<n; i++)
	fprintf(fp, "%le\n", U(points[i][0],points[i][1]));
	fprintf(fp, "\n\n");
	for(int i=0; i<n; i++)
		fprintf(fp, "%le\n", abs(U(points[i][0],points[i][1])-Xk[i]));
	fprintf(fp, "\n\n");
	   for(int i=0; i<n; i++)
norm+=pow((U(points[i][0],points[i][1])-Xk[i]),2);

	   fprintf(fp, "%le\n",sqrt(norm));
}

int main()
{
	setlocale(LC_ALL, "Russian");
	if (!input())
	{
		cout << "Ќеверные данные" << endl;
		system("pause");
		return 0;
	}
   	genSet();
	output_grid();
	create_SLAU();
	//надо убрать эти выводу в файл. они мен€ бес€т))
	output_matrix();
	New_input_matrix_1();
	Jakoby(1);
	WriteRes();
	_getch();
	return 0;
}
