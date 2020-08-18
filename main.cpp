#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <random>
#include <fstream>

using namespace std;
float** read_file(std::string filename,int rows,int cols);
void cluster(int * cluster, float **Pressurelattice, int **copy, int **tab, int nc, int nl, float Pmax);

int main(int argc, char **argv)
{
    int **tab,**copy,nl,nc,a,b,c,d,count=0,fill,*clus,reach=0;
    float **Pressurelattice,alea,gamma,Rmin,Rmax,Pmin,Pc,Pair,P_env,Peau,dP,nbTotal=0,nbAir,satAir,rl,rc,tmp=0,delta_l=0.001,alpha,g0=9.81,delta_rho,rho_eau=1200.0,rho_air=1.25;
    float dimensionalconstant;
    FILE *out,*net,*env;
    clus=(int*) malloc (4*sizeof(int));

    //Lattice dimensions
    nl = 232;
    nc = 330;
    alpha = 1;
    dimensionalconstant = 0.14/nl;
    //cout<<argv[1];

    tab=(int**)malloc((nl+2)*sizeof(int*));
    for (int i=0;i<(nl+2);i++)
    {
        tab[i]=(int*)malloc((nc+2)*sizeof(int));
    }

    for (int i=0;i<(nl+2);i++)
        for (int j=0;j<(nc+2);j++)
        {
            if(i==0)
                tab[i][j]='F'; //Filter
            else
                tab[i][j]='B';
        }

    for (int i=0;i<(nl+2);i++)
    {
        for (int j=0;j<(nc+2);j++)
        {

        }
    }

    //45 degree
    for (int i=1;i<(nl+1);i++)
        if (i%2!=0)
        {
            for (int j=1;j<(nc+1);j++)
            {
                if (j%2==0)
                {
                    tab[i][j]=0; //water phase
                    nbTotal++;
                }
                else
                    tab[i][j]=1; //rock
            }
        }
        else
        {
            for (int j=1;j<(nc+1);j++)
            {
                if (j%2!=0)
                {
                    tab[i][j]=0;
                    nbTotal++;
                }
                else
                    tab[i][j]=1;
            }
        }

    //Adding air as the initial condition
    for (int j=1;j<(nc+1);j++)
        if (tab[nl][j]==0)
            tab[nl][j]=5;

    copy=(int**)malloc((nl+2)*sizeof(int*));
        for (int i=0;i<(nl+2);i++)
            copy[i]=(int*)malloc((nc+2)*sizeof(int));

    //Pressure lattice
    srand(time(NULL));
    gamma=0.073; //surface tension air/water (units N/m)  20°C
    Rmin=0.000240;
    Rmax=0.000730;
    
    //Pmin=2*gamma/Rmax;
    //Pmax=2*gamma/Rmin;

    //Pressure distribution
    // default_random_engine generator;
    //float std_gauss = strtol(argv[4], NULL, 10);
    //std_gauss = 1/std_gauss;
    //normal_distribution<double> distribution(1000,std_gauss);
    //exponential_distribution<double> distribution(std_gauss);

    int rows = 200;
    int cols = 2;
    float** floats;
    if( !(floats = read_file("/home/monem/Dropbox/code/c++/test/poredist.txt",rows,cols) ) ){return 0;}

    Pmin = 200;
    float Pmax = 485.75;
    //rho_eau = strtol(argv[4], NULL, 10);

    Pressurelattice=(float**)malloc(nl*sizeof(int*));
    for (int i=0;i<(nl-1);i++)
    {
        Pressurelattice[i]=(float*)malloc((nc-1)*sizeof(float));
    }

    double diff, epsilon;
    int idx;
    epsilon = 1;

    for (int i=0;i<(nl-1);i++)
    {
        for (int j=0;j<(nc-1);j++)
        {
            alea=(((float)(rand())/(float)(RAND_MAX))); //random number

            for(int k=0;k<rows;k++){
                diff = fabs(alea-floats[1][k]);
                if(diff <= epsilon)
                {
                    epsilon = diff;
                    idx = k;
                }

            }

            //cout<<alea<<" "<<floats[1][idx]<<endl;
            delta_rho=rho_eau-rho_air;
            //Pc=(float)((float)(alea*(Pmax-Pmin)+Pmin)-(dimensionalconstant*i*delta_rho*g0*sin(alpha*atan(1)*4.0/180)));
            Pc=(floats[0][idx] -(dimensionalconstant*i*delta_rho*g0*sin(alpha*atan(1)*4.0/180)));
            Pressurelattice[i][j]=Pc;
            idx = 0;
            epsilon = 1;
        }
    }

    //open file and write.
    out=fopen("PS_test_60.txt","w");
    //out = fopen(argv[2],"w");
    net=fopen("lattice_delta_width_440_angle60.txt","w");
    //net = fopen(argv[1],"w");
    env=fopen("contour60_test.txt","w");
    //env = fopen(argv[3],"w");

    if (nc%2==0)
    {
        rc=(float)(nc/2);
    }
    else
    {
        rc=(float)((nc-1)/2);
    }

    rl=(float)nl;
    //The total number of sites.
    printf("Total number of sites: %f\n",nbTotal);
    cout<<endl;

    //Major IP L00p
    Pair=0;
    Peau=0;
    P_env=0;
    nbAir=(int)(nc/2); //initial number of air-filled pores (on i = nl)
    int db = 1;

    while(1)
    {
        satAir=(float)(nbAir/nbTotal); // saturation en air
        fprintf(out,"%f %.4f\n",Pair,satAir);
        fprintf(env,"%f %.4f\n",P_env,satAir);
        cluster(clus,Pressurelattice,copy,tab,nc,nl,Pmax);

        a=clus[0]; b=clus[1]; c=clus[2]; d=clus[3];

        tab[a][b]=5;
        nbAir++;
        Pc=Pressurelattice[c][d];
        dP=Pc-(Pair-Peau);
        Pair=Pair+dP;

        if (Pair>tmp)
        {
            P_env=Pair;
            tmp=P_env;
        }
        count=0;

        for (int j=1;j<(nc+1);j++)
            if (tab[1][j]==0)
                count++;
        if (a==1) //we touch at the top of the network
        {
            if (reach==0) //if one has never reached the filter
            {
                fprintf(net,"%.4f",satAir);
                fill=rc-2;
                while(fill)
                {
                    fprintf(net,"%4d",1);
                    fill--;
                }
                fprintf(net,"\n");
                printf("Air saturation at the percolation threshold: %f \n",satAir);
                reach=1;
            }
        }

        if (count == 0) //if all the pores of the bottom are top with air
        {
            for (int i=1;i<nl+1;i++)
            {
                for (int j=1;j<nc-1;j++)
                {
                    if (tab[i][j]==1)
                        continue;
                    else
                        fprintf(net,"%4d",tab[i][j]);
                }
                fprintf(net,"\n");
            }
            fprintf(net,"\n");

            printf("Final air saturation : %f\n",satAir);
            break; // Stop when there is more water on the edge
        }
        db++;

    }

    printf("Final air Number of sites: %f\n",nbAir);
    fclose(out);
    fclose(net);
    fclose(env);

    return 0;

}
void cluster(int *cluster, float **Pressurelattice, int **copy, int **tab, int nc, int nl, float Pmax)
{
    int i,j,l,c,m,n,it[2]={-1,1},ip[2]={-2,-1};
    float Plim=Pmax+1,Pc;

    for (i=0;i<(nl+2);i++)
        for (j=0;j<(nc+2);j++)
        {
            if (tab[i][j]=='F')
                copy[i][j]=1; // water is connected to the outside.
            else
                copy[i][j]=0;
        }

    for (i=0;i<nl+2;i++)
        for (j=0;j<nc+2;j++)
        {
            if (copy[i][j]==1) //if there is water from outside,
            {
                for (n=0;n<2;n++) //search nearby sites.
                {
                    for (m=0;m<2;m++)
                    {
                        l=i+it[n];
                        c=j+it[m];
                        if (l>=0 && l<nl+2 && c>=0 && c<nc+2 && tab[l][c]==0)
                        {
                            copy[i+it[n]][j+it[m]]=1; //if the neighboring site also contains water, it is recorded in copy

                        }
                    }
                }
            }
        }

    for (i=1;i<nl+1;i++)
        for (j=1;j<nc+1;j++)
            if (copy[i][j]==1) // if the site contains water from outside,
                    {
                for (m=0;m<2;m++)
                {
                    for (n=0;n<2;n++)
                    {
                        if (tab[i+it[m]][j+it[n]]==5) // if one seeks the neighboring site contains air
                        {
                            Pc=Pressurelattice[i+ip[m]][j+ip[n]]; // Looking at the value of the corresponding threshold pressure Pc
                            if (Pc<Plim)
                            {
                                Plim=Pc; // Search mini (Pc) (to invade the site)
                                cluster[0]=i;
                                cluster[1]=j;
                                cluster[2]=i+ip[m];
                                cluster[3]=j+ip[n];

                            }
                        }
                    }
                }
            }

    return;
}

float** read_file(std::string filename,int rows,int cols)
{
    std::fstream file;
    file.open(filename.c_str(), std::ios::in);
    if(!file.is_open()){return 0;}
    float** floats = new float*[cols+1];
    for(int i = 0; i <cols;++i){ floats[i] = new float[rows+1]; }
    //read each row
    for(int i =0;i<rows;++i)
    {
        for(int j =0;j<cols;++j)//push into the col
        { file >>floats[j][i]; }
    }
    file.close();

    return floats;
}
