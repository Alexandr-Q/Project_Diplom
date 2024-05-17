#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <map>
#include "time.h"

void Num (std::string & Str, std::vector <int> & Vec)           //вспомогательная для чтения из файла
{
const char &rasd = ' '  ;
   Vec.clear();
    std::string TempS = "";
    TempS.push_back(rasd);
    TempS.push_back(rasd);
   int q = 0;
    while (Str.find (TempS, q) != -1)
    {
        q = Str.find (TempS, q);
        Str.erase(q, 1);
    }
    while (Str[0] == rasd)
    {Str.erase(0, 1);}
    while (Str[Str.length()-1] == rasd)
    {Str.erase((Str.length()-1), 1);}
    std::string TS = TempS;
    TS.pop_back();
    TempS = "";
    int b=0;
    int e=0;
    int r;
    while (Str.find (TS, b) != -1)
    {
        e = Str.find (TS, b)-1;
        TempS = Str.substr(b, e-b+1);
         r = atoi(TempS.c_str());
         Vec.push_back(r);
        b = e+2;
        TempS.clear();
    }
    TempS = Str.substr(b, Str.length()-b);
    r = atoi(TempS.c_str());
    Vec.push_back(r);
    TempS.clear();
}

void NumMet (std::string & Str, int & Nomer, std::string & Metka)           //вспомогательная для чтения из файла меток
{
const char &rasd = ' '  ;
    std::string TempS = "";
    TempS.push_back(rasd);
    TempS.push_back(rasd);
   int q = 0;
    while (Str.find (TempS, q) != -1)
    {
        q = Str.find (TempS, q);
        Str.erase(q, 1);
    }
    while (Str[0] == rasd)
    {Str.erase(0, 1);}
    while (Str[Str.length()-1] == rasd)
    {Str.erase((Str.length()-1), 1);}
    std::string TS = TempS;
    TS.pop_back();
    TempS = "";
    int b=0;
    int e=0;
    int r;
    while (Str.find (TS, b) != -1)
    {
        e = Str.find (TS, b)-1;
        TempS = Str.substr(b, e-b+1);
         r = atoi(TempS.c_str());
         Nomer = r;
        b = e+2;
        TempS.clear();
    }
    TempS = Str.substr(b, Str.length()-b);
    Metka = TempS;
    TempS.clear();
}

int FileToVec (std::ifstream & FileData, std::vector <int> & VecData)       // чтение графа из файла
{
    std::string TempData = "";
    VecData.clear();
    std::vector <int> TempVec;
    while (!FileData.eof())
    {
       TempVec.clear();
        getline (FileData, TempData);
        if (TempData.length()!=0)
        {
             Num(TempData, TempVec);
             if (TempVec.size()!=2) {VecData.clear(); return -1;}
             VecData.push_back(TempVec[0]);
             VecData.push_back(TempVec[1]);
        }
    }
    if (VecData.size()==0) return -1;
    return 0;
}

int FileToMet (std::ifstream & FileData, std::map<int, std::string> & MetData)  // чтение меток из файла
{
    std::string TempData = "";
    MetData.clear();
    std::string TempStr;
    int TempInt;
    while (!FileData.eof())
    {
       TempStr.clear();
       TempInt = 0;
        getline (FileData, TempData);
        if (TempData.length()!=0)
        {
             NumMet(TempData, TempInt, TempStr);
             MetData.emplace(std::make_pair(TempInt, TempStr));
        }
    }
    if (MetData.size()==0) return -1;
    return 0;
}

int VecToFile (const std::vector <int> &VecRes, std::ofstream &FileRes)     // запись результата в файл
{
    for (int i=0; i<VecRes.size();i+=2)
    {
        FileRes<< VecRes[i]<<' '<<VecRes[i+1];
         FileRes<< std::endl;
    }
    FileRes<< std::endl;
    return 0;
}

int PerNum (std::vector <int> & Vec)        // вспомогательная: перенумерация вершин графа
{
    std::set<int> AA;
    AA.clear();
    std::map <int, int> P;
    P.clear();
    int s=0;
    for (int i=0; i<Vec.size(); i++)
    {
        AA.insert(Vec[i]);
    }
    for (auto i=AA.begin(); i!=AA.end(); i++)
    {
        P.insert(std::pair<int, int>(*i, s));
        s++;
    }
    for (int i=0; i<Vec.size(); i++)
    {
        Vec[i] = P[Vec[i]];
    }
    return 0;
}

int KratReber (const std::vector <int> &Graph, std::map <std::pair < int, int>, int> &Krat, bool npr)   // вспомогательная: подсчет кратности ребер
{
    Krat.clear();
    std::pair < int, int> R, RR ;
    std::pair < std::pair < int, int>, int> K, KK;
    if (npr)
    {
        for (int i=0; i<Graph.size(); i+=2)
        {
            R = std::make_pair(Graph[i], Graph[i+1]);
            K = std::make_pair(R, 1);
            if (Krat.find(R)!=Krat.end())
            {
                Krat[R] = Krat[R]+1;
                continue;
            }
            if (Krat.find(R)==Krat.end())
            {
                Krat.insert(K);
                continue;
            }
        }
    }
    else
    {
        for (int i=0; i<Graph.size(); i+=2)
        {
            R = std::make_pair(Graph[i], Graph[i+1]);
            K = std::make_pair(R, 1);
            RR = std::make_pair(Graph[i+1], Graph[i]);
            KK = std::make_pair(RR, 1);
            if (Krat.find(R)!=Krat.end())
            {
                Krat[R] = Krat[R]+1;
                Krat[RR] = Krat[RR]+1;
                continue;
            }
            if (Krat.find(R)==Krat.end())
            {
                Krat.insert(K);
                Krat.insert(KK);
                continue;
            }
        }
    }
    return 0;
}

int RasVer (std::vector <int> &Graph, std::vector <int> & Rass, int v, int M = INT_MIN) // вспомогательная: подсчет расстояний от заданной вершины
{
    Rass.clear();
    int nmin, nmax, flag;
    nmin = INT_MAX;
    nmax = INT_MIN;
    if ((M<0)&&(M != INT_MIN)) return -1;
    for (int ii=0; ii<Graph.size(); ii++)
    {
        if (Graph[ii]>nmax) nmax = Graph[ii];
        if (Graph[ii]<nmin) nmin = Graph[ii];
    }
    if (nmax>M) M = nmax;
    Rass.resize(M+1, INT_MAX);
    Rass[v] = 0;
    flag=0;
    for (int s=nmin; s<=nmax; s++)
    {
        flag=0;
        for (int j=0; j<Graph.size(); j=j+2)
        {
            if (Rass[Graph[j]]== INT_MAX) continue;
            if ( Rass[Graph[j+1]] > Rass[Graph[j]] + 1)
            {
                Rass[Graph[j+1]] = Rass[Graph[j]] + 1;
                flag=1;
            }
        }
        if (flag==0)
        {
            break;
        }
        if (s==nmax)
        {
            flag=-1;
            break;
        }
    }
    if (flag==-1)
    {
        Rass.clear();
        return -1;
    }
    return 0;
}

int ZamVec (std::vector <int> & V, int v1, int v2)    // вспомогательная: замена элементов вектора
{
    int t = V[v1];
    V[v1] = V[v2];
    V[v2] = t;
    return 0;
}

double TimeI = 0;

int PoiskGraph (std::vector <int> A, std::vector <int> B, std::set<std::vector <int>> & A_B, bool npr = true,
                std::map<int, std::string> mgA= {}, std::map<int, std::string> mgB= {})                         // основная функция
{
    bool NM=true;                                   //наличие меток
    bool QM = false;                                //корректность меток
    int sB=B.size();                                // количество вершин В
    int a, f1, f2, f3, dl, ll, sx, sz,  t1, t2 =0;  //вспомогательные переменные
    int minA, minB  = INT_MAX;                      //номера мин. вершин
    int maxA, maxB = INT_MIN;                       //номера макс. вершин

    A_B.clear();                                    //найденные подграфы
    if (A.size()==0) return -1;
    if (B.size()==0) return -1;
    if (B.size()>A.size()) return -1;
    if ((A.size() % 2)==1) return -1;
    if ((B.size() % 2)==1) return -1;
    if (mgA.size()<mgB.size()) return -1;

    if ((mgA.size()==0)&&(mgB.size()==0))
        NM=false;

    if (mgB.size()==0)
    {
        PerNum (B);
    }

    if (mgA.size()>0)
        for (auto iA=mgA.begin(); iA!=mgA.end(); iA++)
        {
            for (int i = 0; i<A.size(); i++)
            {
                if (A[i] == iA->first) QM = true;
            }
            if (QM==false) return -1;
            QM = false;
        }

    if (mgB.size()>0)
        for (auto iB=mgB.begin(); iB!=mgB.end(); iB++)
        {
            for (int i = 0; i<B.size(); i++)
            {
                if (B[i] == iB->first) QM = true;
            }
            if (QM==false) return -1;
            QM = false;
        }

        double Time = time(NULL);




    for (int ii=0; ii<A.size(); ii++)
    {
        if (A[ii]>maxA) maxA = A[ii];
        if (A[ii]<minA) minA = A[ii];
    }

    if (minA<0) return -1;

    for (int ii=0; ii<B.size(); ii++)
    {
        if (B[ii]>maxB) maxB = B[ii];
        if (B[ii]<minB) minB = B[ii];
    }

    if (minB<0) return -1;

    std::vector <int> AA;
    AA.clear();
    std::vector <int> BB;
    BB.clear();
    std::vector <int> GR;
    GR.clear();  //
    std::vector <int> MS;
    MS.clear();
    std::map <std::pair < int, int>, int> KratA;
    std::map <std::pair < int, int>, int> KratB;
    KratA.clear();
    KratB.clear();

    if (!npr)
    {
        for (int i=0; i<A.size(); i+=2)
        {
            if (A[i]>A[i+1]) ZamVec(A, i, i+1);
        }
        for (int i=0; i<B.size(); i+=2)
        {
            if (B[i]>B[i+1]) ZamVec(B, i, i+1);
        }
    }

    KratReber (A, KratA, npr);
    KratReber (B, KratB, npr);
    A.clear();
    B.clear();

    for (auto it = KratA.begin(); it!=KratA.end(); it++)
    {
        {
            A.push_back((it->first).first);
            A.push_back((it->first).second);
        }
    }

    for (auto it = KratB.begin(); it!=KratB.end(); it++)
    {
        {
            B.push_back((it->first).first);
            B.push_back((it->first).second);
        }
    }

    std::vector <int> PSZA;                     // полустепени захода вершин графа А
    std::vector <int> PSVA;                     // полустепени исхода вершин графа А
    std::vector <int> PSZB;                     // для графа В
    std::vector <int> PSVB;                     // для графа В
    std::vector <std::vector<int> >RofA;        // ребра графа А
    std::vector <std::vector<int> >RofB;        // ребра графа В
    std::vector <std::vector<int>> RofBfromA (RofB.size()); // соответствие ребер А и В
    std::vector <int> D;                        // вспомогательный для подсчета расстояний между вершинами
    std::vector <std::vector<int>> DA;          // кратчайшее расстояние между вершинами А
    std::vector <std::vector <int>> DB;         // кратчайшее расстояние между вершинами В
    std::vector <int> E;                        // вспомогательный для дублирования ребер ненаправленного графа
    DA.clear();
    DB.clear();
    D.clear();
    E.clear();

    minA = INT_MAX;
    maxA = INT_MIN;

    for (int ii=0; ii<A.size(); ii++)
    {
        if (A[ii]>maxA) maxA = A[ii];
        if (A[ii]<minA) minA = A[ii];
    }

    PSZA.clear();
    PSZA.resize(maxA+1, 0);
    PSVA.clear();
    PSVA.resize(maxA+1, 0);


    RofA.clear();                               // для графа А

    for (int a = 0; a<A.size(); a+=2)
    {
        RofA.push_back({A[a], A[a+1]});
        if (!npr)
            RofA.push_back({A[a+1],A[a]});
    }
    for (int a = 0; a<A.size(); a+=2)
    {
        if (npr)
        {
            PSZA[(A[a+1])]++;
            PSVA[(A[a])]++;
        }
        if (!npr)
        {
            PSZA[(A[a+1])]++;
            PSZA[(A[a])]++;
        }
    }

    RofB.clear();                               // для графа В

    for (int a = 0; a<B.size(); a+=2)
    {
        RofB.push_back({B[a], B[a+1]});
    }

    PSZB.clear();
    PSZB.resize(maxB+1, 0);
    PSVB.clear();
    PSVB.resize(maxB+1, 0);

    for (int a = 0; a<B.size(); a+=2)
    {
        if (npr)
        {
            PSZB[(B[a+1])]++;
            PSVB[(B[a])]++;
        }
        if (!npr)
        {
            PSZB[(B[a+1])]++;
            PSZB[(B[a])]++;
        }
    }
    sz= RofA.size();

    RofBfromA.clear();
    RofBfromA.resize(RofB.size());              // сравнение А и В

    for (int j=0; j<RofB.size(); j++)
        RofBfromA[j].clear();

    t1 = RofA.size();

    for (int j=0; j<RofB.size(); j++)           // сравнение внутренних
    {
        if ( (PSZB[(RofB[j][0])]+PSVB[RofB[j][0]])==1)
            continue;
        if ( (PSZB[(RofB[j][RofB[j].size()-1])]+PSVB[RofB[j][RofB[j].size()-1]])==1)
            continue;
        if ( (RofB[j][0] == RofB[j][RofB[j].size()-1]) && ((PSZB[(RofB[j][0])]+PSVB[RofB[j][0]])==2) )
            continue;
        for (int i=0; i<RofA.size(); i++)
        {
            if (  ( (RofB[j][0] == RofB[j][RofB[j].size()-1]) && (RofA[i][0] == RofA[i][RofA[i].size()-1]) )
                || ( (RofB[j][0] != RofB[j][RofB[j].size()-1]) && (RofA[i][0] != RofA[i][RofA[i].size()-1]) ))
                if ( (PSZB[(RofB[j][0])]<= PSZA[(RofA[i][0])]) && (PSZB[(RofB[j][RofB[j].size()-1])]<= PSZA[(RofA[i][RofA[i].size()-1])])
                    && (PSVB[(RofB[j][0])]<= PSVA[(RofA[i][0])]) && (PSVB[(RofB[j][RofB[j].size()-1])]<= PSVA[(RofA[i][RofA[i].size()-1])]))
                {
                    f2=0;
                    for (int a=1; a<RofB[j].size(); a++)
                    {
                       if (npr)
                            if (KratA[std::pair <int, int>(RofA[i][a-1], RofA[i][a]) ] < KratB[std::pair <int, int>(RofB[j][a-1], RofB[j][a])])
                            {
                                f2=1;
                                break;
                            }
                        if (!npr)
                            if (KratA[std::pair <int, int>( std::min(RofA[i][a-1], RofA[i][a]), std::max(RofA[i][a-1], RofA[i][a])) ]
                                < KratB[std::pair <int, int>( std::min(RofB[j][a-1], RofB[j][a]), std::max(RofB[j][a-1], RofB[j][a]))])
                            {
                                f2=1;
                                break;
                            }
                        if (NM==true)                           // проверка меток
                        {
                            if (mgB.find(RofB[j][a-1])!=mgB.end())
                            {
                                if (mgA.find(RofA[i][a-1])!=mgA.end())
                                    if (mgB[RofB[j][a-1]]!=mgA[RofA[i][a-1]])
                                    {
                                        f2=1;
                                        break;
                                    }
                                if (mgA.find(RofA[i][a-1])==mgA.end())
                                {
                                    f2=1;
                                    break;
                                }
                            }
                            if (mgB.find(RofB[j][a-1])==mgB.end())
                            {
                                if (mgA.find(RofA[i][a-1])!=mgA.end())
                                {
                                    f2=1;
                                    break;
                                }
                            }
                            if (mgB.find(RofB[j][a])!=mgB.end())
                            {
                                if (mgA.find(RofA[i][a])!=mgA.end())
                                    if (mgB[RofB[j][a]]!=mgA[RofA[i][a]])
                                    {
                                        f2=1;
                                        break;
                                    }
                                if (mgA.find(RofA[i][a])==mgA.end())
                                {
                                    f2=1;
                                    break;
                                }
                            }
                            if (mgB.find(RofB[j][a])==mgB.end())
                            {
                                if (mgA.find(RofA[i][a])!=mgA.end())
                                {
                                    f2=1;
                                    break;
                                }
                            }
                        }

                    }
                    if (f2==0)
                    {
                        RofBfromA[j].push_back(i);
                        RofBfromA[j].push_back(0);
                        RofBfromA[j].push_back(RofB[j].size());
                    }
                }
        }
    }

    for (int j=0; j<RofB.size(); j++)           //сравнение концевых
    {
        if( ( (PSZB[(RofB[j][0])]+PSVB[RofB[j][0]])==1) && ((PSZB[(RofB[j][RofB[j].size()-1])]+PSVB[RofB[j][RofB[j].size()-1]])==1))
            continue;
        if(( ( (PSZB[(RofB[j][0])]+PSVB[RofB[j][0]])==1)
                && ((PSZB[(RofB[j][RofB[j].size()-1])]+PSVB[RofB[j][RofB[j].size()-1]])>1)) || ( ( (PSZB[(RofB[j][0])]+PSVB[RofB[j][0]])>1)
                        && ((PSZB[(RofB[j][RofB[j].size()-1])]+PSVB[RofB[j][RofB[j].size()-1]])==1)))
            for (int i=0; i<RofA.size(); i++)
            {
                if ( (PSZB[(RofB[j][0])]<= PSZA[(RofA[i][0])]) && (PSZB[(RofB[j][RofB[j].size()-1])]<= PSZA[(RofA[i][RofA[i].size()-1])])
                        && (PSVB[(RofB[j][0])]<= PSVA[(RofA[i][0])]) && (PSVB[(RofB[j][RofB[j].size()-1])]<= PSVA[(RofA[i][RofA[i].size()-1])]))
                {
                    if ((PSZB[(RofB[j][RofB[j].size()-1])]+PSVB[RofB[j][RofB[j].size()-1]])==1)
                    {
                        f2=0;
                        for (int a=1; a<RofB[j].size(); a++)
                        {
                            if (npr)
                                if (KratA[std::pair <int, int>(RofA[i][a-1], RofA[i][a]) ] < KratB[std::pair <int, int>(RofB[j][a-1], RofB[j][a])])
                                {
                                    f2=1;
                                    break;
                                }
                            if (!npr)
                                if (KratA[std::pair <int, int>( std::min(RofA[i][a-1], RofA[i][a]), std::max(RofA[i][a-1], RofA[i][a])) ]
                                        < KratB[std::pair <int, int>( std::min(RofB[j][a-1], RofB[j][a]), std::max(RofB[j][a-1], RofB[j][a]))])
                                {
                                    f2=1;
                                    break;
                                }
                            if (NM==true)
                            {
                                if (mgB.find(RofB[j][a-1])!=mgB.end())
                                {
                                    if (mgA.find(RofA[i][a-1])!=mgA.end())
                                        if (mgB[RofB[j][a-1]]!=mgA[RofA[i][a-1]])
                                        {
                                            f2=1;
                                            break;
                                        }
                                    if (mgA.find(RofA[i][a-1])==mgA.end())
                                    {
                                        f2=1;
                                        break;
                                    }
                                }
                                if (mgB.find(RofB[j][a-1])==mgB.end())
                                {
                                    if (mgA.find(RofA[i][a-1])!=mgA.end())
                                    {
                                        f2=1;
                                        break;
                                    }
                                }
                                if (mgB.find(RofB[j][a])!=mgB.end())
                                {
                                    if (mgA.find(RofA[i][a])!=mgA.end())
                                        if (mgB[RofB[j][a]]!=mgA[RofA[i][a]])
                                        {
                                            f2=1;
                                            break;
                                        }
                                    if (mgA.find(RofA[i][a])==mgA.end())
                                    {
                                        f2=1;
                                        break;
                                    }
                                }
                                if (mgB.find(RofB[j][a])==mgB.end())
                                {
                                    if (mgA.find(RofA[i][a])!=mgA.end())
                                    {
                                        f2=1;
                                        break;
                                    }
                                }
                            }
                        }
                        if (f2==0)
                        {
                            RofBfromA[j].push_back(i);
                            RofBfromA[j].push_back(0);
                            RofBfromA[j].push_back(RofB[j].size());
                        }
                    }
                    if ( (PSZB[(RofB[j][0])]+PSVB[RofB[j][0]])==1)
                    {
                        f2=0;
                        for (int a=1; a<RofB[j].size(); a++)
                            {
                            if (npr)
                                if (KratA[std::pair <int, int>(RofA[i][RofA[i].size()-RofB[j].size()+a-1], RofA[i][RofA[i].size()-RofB[j].size()+a])]
                                        < KratB[std::pair <int, int>(RofB[j][a-1], RofB[j][a])])
                                {
                                    f2=1;
                                    break;
                                }
                            if (!npr)
                                if (KratA[std::pair <int, int>( std::min(RofA[i][RofA[i].size()-RofB[j].size()+a-1],RofA[i][RofA[i].size()-RofB[j].size()+a]),
                                                                std::max(RofA[i][RofA[i].size()-RofB[j].size()+a-1], RofA[i][RofA[i].size()-RofB[j].size()+a]))]
                                        < KratB[std::pair <int, int>( std::min(RofB[j][a-1], RofB[j][a]), std::max(RofB[j][a-1], RofB[j][a]))])
                                {
                                    f2=1;
                                    break;
                                }
                            if (NM==true)
                            {
                                if (mgB.find(RofB[j][a-1])!=mgB.end())
                                {
                                    if (mgA.find(RofA[i][RofA[i].size()-RofB[j].size()+a-1])!=mgA.end())
                                        if (mgB[RofB[j][a-1]]!=mgA[RofA[i][RofA[i].size()-RofB[j].size()+a-1]])
                                        {
                                            f2=1;
                                            break;
                                        }
                                    if (mgA.find(RofA[i][RofA[i].size()-RofB[j].size()+a-1])==mgA.end())
                                    {
                                        f2=1;
                                        break;
                                    }
                                }
                                if (mgB.find(RofB[j][a-1])==mgB.end())
                                {
                                    if (mgA.find(RofA[i][RofA[i].size()-RofB[j].size()+a-1])!=mgA.end())
                                    {
                                        f2=1;
                                        break;
                                    }
                                }
                                if (mgB.find(RofB[j][a])!=mgB.end())
                                {
                                    if (mgA.find(RofA[i][RofA[i].size()-RofB[j].size()+a])!=mgA.end())
                                        if (mgB[RofB[j][a]]!=mgA[RofA[i][RofA[i].size()-RofB[j].size()+a]])
                                        {
                                            f2=1;
                                            break;
                                        }
                                    if (mgA.find(RofA[i][RofA[i].size()-RofB[j].size()+a])==mgA.end())
                                    {
                                        f2=1;
                                        break;
                                    }
                                }
                                if (mgB.find(RofB[j][a])==mgB.end())
                                {
                                    if (mgA.find(RofA[i][RofA[i].size()-RofB[j].size()+a])!=mgA.end())
                                    {
                                        f2=1;
                                        break;
                                    }
                                }
                            }
                        }
                        if (f2==0)
                        {
                            RofBfromA[j].push_back(i);
                            RofBfromA[j].push_back(RofA[i].size()-RofB[j].size());
                            RofBfromA[j].push_back(RofB[j].size());
                        }
                    }
                }
            }
    }

    for (int j=0; j<RofB.size(); j++)       //сравнение изолированных
    {
        if( ( (PSZB[(RofB[j][0])]+PSVB[RofB[j][0]])==1) && ((PSZB[(RofB[j][RofB[j].size()-1])]+PSVB[RofB[j][RofB[j].size()-1]])==1))
            for (int i=0; i<t1; i++)
            {
                if ( (RofA[i][0] != RofA[i][RofA[i].size()-1]))
                {
                    for (int z=0; z<=RofA[i].size()-RofB[j].size(); z++)
                    {
                        f2=0;
                        for (int a=1; a<RofB[j].size(); a++)
                            {
                            if (npr)
                                if (KratA[std::pair <int, int>(RofA[i][z+a-1], RofA[i][z+a]) ] < KratB[std::pair <int, int>(RofB[j][a-1], RofB[j][a])])
                                {
                                    f2=1;
                                    break;
                                }
                            if (!npr)
                                if (KratA[std::pair <int, int>( std::min(RofA[i][z+a-1], RofA[i][z+a]), std::max(RofA[i][z+a-1], RofA[i][z+a]))]
                                        < KratB[std::pair <int, int>( std::min(RofB[j][a-1], RofB[j][a]), std::max(RofB[j][a-1], RofB[j][a]))])
                                {
                                    f2=1;
                                    break;
                                }
                            if (NM==true)
                            {
                                if (mgB.find(RofB[j][a-1])!=mgB.end())
                                {
                                    if (mgA.find(RofA[i][z+a-1])!=mgA.end())
                                        if (mgB[RofB[j][a-1]]!=mgA[RofA[i][z+a-1]])
                                        {
                                            f2=1;
                                            break;
                                        }
                                    if (mgA.find(RofA[i][z+a-1])==mgA.end())
                                    {
                                        f2=1;
                                        break;
                                    }
                                }
                                if (mgB.find(RofB[j][a-1])==mgB.end())
                                {
                                    if (mgA.find(RofA[i][z+a-1])!=mgA.end())
                                    {
                                        f2=1;
                                        break;
                                    }
                                }
                                if (mgB.find(RofB[j][a])!=mgB.end())
                                {
                                    if (mgA.find(RofA[i][z+a])!=mgA.end())
                                        if (mgB[RofB[j][a]]!=mgA[RofA[i][z+a]])
                                        {
                                            f2=1;
                                            break;
                                        }
                                    if (mgA.find(RofA[i][z+a])==mgA.end())
                                    {
                                        f2=1;
                                        break;
                                    }
                                }
                                if (mgB.find(RofB[j][a])==mgB.end())
                                {
                                    if (mgA.find(RofA[i][z+a])!=mgA.end())
                                    {
                                        f2=1;
                                        break;
                                    }
                                }
                            }
                        }
                        if (f2==0)
                        {
                            RofBfromA[j].push_back(i);
                            RofBfromA[j].push_back(z);
                            RofBfromA[j].push_back(RofB[j].size());
                        }
                    }
                }
            }
    }

    for (int j=0; j<RofB.size(); j++)                       //сравнение циклов
    {
        if( (RofB[j][0] != RofB[j][RofB[j].size()-1]))
            continue;
        for (int i=0; i<t1; i++)
        {
            if ( (RofA[i][0] == RofA[i][RofA[i].size()-1]))
            {
                f2=0;
                for (int a=1; a<RofB[j].size(); a++)
                    {
                    if (npr)
                        if (KratA[std::pair <int, int>(RofA[i][a-1], RofA[i][a]) ] < KratB[std::pair <int, int>(RofB[j][a-1], RofB[j][a])])
                        {
                            f2=1;
                            break;
                        }
                    if (!npr)
                        if (KratA[std::pair <int, int>( std::min(RofA[i][a-1], RofA[i][a]), std::max(RofA[i][a-1], RofA[i][a])) ]
                            < KratB[std::pair <int, int>( std::min(RofB[j][a-1], RofB[j][a]), std::max(RofB[j][a-1], RofB[j][a]))])
                        {
                            f2=1;
                            break;
                        }
                    if (NM==true)
                    {
                        if (mgB.find(RofB[j][a-1])!=mgB.end())
                        {
                            if (mgA.find(RofA[i][a-1])!=mgA.end())
                                if (mgB[RofB[j][a-1]]!=mgA[RofA[i][a-1]])
                                {
                                    f2=1;
                                    break;
                                }
                            if (mgA.find(RofA[i][a-1])==mgA.end())
                            {
                                f2=1;
                                break;
                            }
                        }
                        if (mgB.find(RofB[j][a-1])==mgB.end())
                        {
                            if (mgA.find(RofA[i][a-1])!=mgA.end())
                            {
                                f2=1;
                                break;
                            }
                        }
                        if (mgB.find(RofB[j][a])!=mgB.end())
                        {
                            if (mgA.find(RofA[i][a])!=mgA.end())
                                if (mgB[RofB[j][a]]!=mgA[RofA[i][a]])
                                {
                                    f2=1;
                                    break;
                                }
                            if (mgA.find(RofA[i][a])==mgA.end())
                            {
                                f2=1;
                                break;
                            }
                        }
                        if (mgB.find(RofB[j][a])==mgB.end())
                        {
                            if (mgA.find(RofA[i][a])!=mgA.end())
                            {
                                f2=1;
                                break;
                            }
                        }
                    }
                }
                if (f2==0)
                {
                    RofBfromA[j].push_back(i);
                    RofBfromA[j].push_back(0);
                    RofBfromA[j].push_back(RofB[j].size());
                }
            }
        }
    }

    for (int i=0; i<RofBfromA.size(); i++)
    {
        if (RofBfromA[i].size()==0)
        {
            return 0;                   //если нет соответствий
        }
    }

    D.clear();

    if ((DA.size()==0)||(DB.size()==0))
    {
        sz=A.size();
        if (!npr)
        {
            E=A;
            std::reverse(E.begin(), E.end());
            A.reserve(A.size()+E.size());
            A.insert(A.end(), E.begin(), E.end());
            E.clear();
        }
        for (int i=0; i<=maxA; i++)     //подсчет мин. расстояний между вершинами А
        {
            RasVer (A, D, i);
            DA.push_back(D);
        }
        sz=B.size();
        if (!npr)
        {
            E=B;
            std::reverse(E.begin(), E.end());
            B.reserve(B.size()+E.size());
            B.insert(B.end(), E.begin(), E.end());
            E.clear();
        }
        for (int i=0; i<=maxB; i++)     //подсчет мин. расстояний между вершинами В
        {
            RasVer (B, D, i);
            DB.push_back(D);
        }
        sz=A.size()/2;
        if (!npr)
            A.resize(sz);
        sz=B.size()/2;
        if (!npr)
            B.resize(sz);
    }

    sz=0;
    sx=INT_MAX;
        for (int i=0; i<RofBfromA.size(); i++)
            if (RofBfromA[i].size()<sx)
            {
                sx= RofBfromA[i].size();
                t2=i;
            }
        f2=0;
        for (int i=0; i<RofBfromA.size(); i++)    //сравнение и сокращение по мин. расстояниям
        {
            if (i==t2) continue;
            for (int y=0; y<RofBfromA[i].size(); y+=3)
            {
                f1=1;
                for (int g=0; g<RofBfromA[t2].size(); g+=3)
                {
                    if ( ((  (DA[RofA[RofBfromA[t2][g]][RofBfromA[t2][g+1]]][RofA[RofBfromA[i][y]][RofBfromA[i][y+1]]])
                             <=  DB[RofB[t2][0]][RofB[i][0]])
                            && (  (DA[RofA[RofBfromA[i][y]][RofBfromA[i][y+1]]][RofA[RofBfromA[t2][g]][RofBfromA[t2][g+1]]])
                                  <=  DB[RofB[i][0]][RofB[t2][0]]) )
                            &&   ((  (DA[RofA[RofBfromA[t2][g]][RofBfromA[t2][g+1]+RofB[t2].size()-1]][RofA[RofBfromA[i][y]][RofBfromA[i][y+1]+RofB[i].size()-1]])
                                     <=  DB[RofB[t2][RofB[t2].size()-1]][RofB[i][RofB[i].size()-1]])
                                  && (  (DA[RofA[RofBfromA[i][y]][RofBfromA[i][y+1]+RofB[i].size()-1]][RofA[RofBfromA[t2][g]][RofBfromA[t2][g+1]+RofB[t2].size()-1]])
                                        <=  DB[RofB[i][RofB[i].size()-1]][RofB[t2][RofB[t2].size()-1]])))
                    {
                        f1=0;
                        break;
                    }
                }
                if (f1==1)
                {
                    f2=1;
                    if (y!=RofBfromA[i].size()-3)
                    {
                        ZamVec(RofBfromA[i], y, RofBfromA[i].size()-3);
                        ZamVec(RofBfromA[i], y+1, RofBfromA[i].size()-2);
                        ZamVec(RofBfromA[i], y+2, RofBfromA[i].size()-1);
                    }
                    RofBfromA[i].pop_back();
                    RofBfromA[i].pop_back();
                    RofBfromA[i].pop_back();
                    y-=3;
                }
            }
        }
        for (int tt=0; tt<RofBfromA.size(); tt++)
        {
            if (RofBfromA[tt].size()==0)
            {
                sz=-1;
                break;
            }
            for (int i=0; i<RofBfromA.size(); i++)
            {
                if (RofBfromA[i].size()==0)
                {
                    sz=-1;
                    break;
                }
                if (i==tt) continue;
                for (int y=0; y<RofBfromA[i].size(); y+=3)
                {
                    f1=1;
                    for (int g=0; g<RofBfromA[tt].size(); g+=3)
                    {
                        if ( ((  (DA[RofA[RofBfromA[tt][g]][RofBfromA[tt][g+1]]][RofA[RofBfromA[i][y]][RofBfromA[i][y+1]]])
                                 <=  DB[RofB[tt][0]][RofB[i][0]])
                                && (  (DA[RofA[RofBfromA[i][y]][RofBfromA[i][y+1]]][RofA[RofBfromA[tt][g]][RofBfromA[tt][g+1]]])
                                      <=  DB[RofB[i][0]][RofB[tt][0]]) )
                                &&   ((  (DA[RofA[RofBfromA[tt][g]][RofBfromA[tt][g+1]+RofB[tt].size()-1]][RofA[RofBfromA[i][y]][RofBfromA[i][y+1]+RofB[i].size()-1]])
                                         <=  DB[RofB[tt][RofB[tt].size()-1]][RofB[i][RofB[i].size()-1]])
                                      && (  (DA[RofA[RofBfromA[i][y]][RofBfromA[i][y+1]+RofB[i].size()-1]][RofA[RofBfromA[tt][g]][RofBfromA[tt][g+1]+RofB[tt].size()-1]])
                                            <=  DB[RofB[i][RofB[i].size()-1]][RofB[tt][RofB[tt].size()-1]])))
                        {
                            f1=0;
                            break;
                        }
                    }
                    if (f1==1)
                    {
                        f2=1;
                        if (y!=RofBfromA[i].size()-3)
                        {
                            ZamVec(RofBfromA[i], y, RofBfromA[i].size()-3);
                            ZamVec(RofBfromA[i], y+1, RofBfromA[i].size()-2);
                            ZamVec(RofBfromA[i], y+2, RofBfromA[i].size()-1);
                        }
                        RofBfromA[i].pop_back();
                        RofBfromA[i].pop_back();
                        RofBfromA[i].pop_back();
                        y-=3;
                    }
                }
                if (RofBfromA[i].size()==0)
                {
                    sz=-1;
                    break;
                }
            }
            if (sz==-1) break;
        }
        for (int i=0; i<RofBfromA.size(); i++)
        {
            if (RofBfromA[i].size()==0)
            {
                sz=-1;
                break;
            }
            for (int y=0; y<RofBfromA[i].size(); y+=3)
            {

                if (  ( (RofB[i][0]==RofB[i][RofB[i].size()-1])
                        && (RofA[RofBfromA[i][y]][RofBfromA[i][y+1]]!= RofA[RofBfromA[i][y]][RofBfromA[i][y+1]-1+RofBfromA[i][y+2]]))
                        || ( (RofB[i][0]!=RofB[i][RofB[i].size()-1])
                             && (RofA[RofBfromA[i][y]][RofBfromA[i][y+1]]== RofA[RofBfromA[i][y]][RofBfromA[i][y+1]-1+RofBfromA[i][y+2]])))
                {
                    f2=1;
                    if (y!=RofBfromA[i].size()-3)
                    {
                        ZamVec(RofBfromA[i], y, RofBfromA[i].size()-3);
                        ZamVec(RofBfromA[i], y+1, RofBfromA[i].size()-2);
                        ZamVec(RofBfromA[i], y+2, RofBfromA[i].size()-1);
                    }
                    RofBfromA[i].pop_back();
                    RofBfromA[i].pop_back();
                    RofBfromA[i].pop_back();
                    y-=3;
                }
                if (RofBfromA[i].size()==0)
                {
                    sz=-1;
                    break;
                }
            }
            if (RofBfromA[i].size()==0)
            {
                sz=-1;
                break;
            }
        }

    sx=0;
    dl= RofB.size();

    std::vector <int> VecL (RofB.size());       // финальное сравнение
    std::vector <int> VecM (RofB.size());
    std::vector<std::vector <int>> MatrB;
    std::vector <int> NumB;

    NumB.clear();
    NumB.resize(maxB+1, -1);

    for (int s=0; s<RofB.size(); s++)
    {
        VecL [s] = 0;
        VecM [s] = RofBfromA[s].size()/3-1;
    }
    VecL[VecL.size()-1]--;

    MatrB.resize (maxB+1-minB);

    for (int row = 0; row< (maxB+1-minB); row++)
    {
        MatrB [row].resize(maxB+1-minB);
        for (int column = 0; column < (maxB+1-minB); column++)
        {
            MatrB [row] [column] = 0;
        }
    }

    for (int q1=0; q1<RofB.size(); q1++)
    {
        if (NumB[RofB[q1][0]]==-1)
        {
            MS.push_back(RofB[q1][0]);
            NumB[RofB[q1][0]] = MS.size()-1;
        }
        for (int q2=0; q2<RofB[q1].size()-1; q2++)
        {
            BB.push_back(RofB[q1][q2]);
            BB.push_back(RofB[q1][q2+1]);
            if (NumB[RofB[q1][q2+1]]==-1)
            {
                MS.push_back(RofB[q1][q2+1]);
                NumB[RofB[q1][q2+1]] = MS.size()-1;
            }
            MatrB  [NumB[RofB[q1][q2]]] [NumB[RofB[q1][q2+1]]]=1;
            if (!npr) MatrB [NumB[RofB[q1][q2+1]]] [NumB[RofB[q1][q2]]] = MatrB  [NumB[RofB[q1][q2]]] [NumB[RofB[q1][q2+1]]];
        }
    }
    std::multiset<std::pair <int, int>> RRR;
    RRR.clear();
    std::set<std::pair <int, int>> Y;
    Y.clear();
    std::vector<int> NumA;
    for (int i = 0; i <= maxA; i++)
        NumA.push_back(-1);
    AA.resize(BB.size());
    GR.resize(B.size());

    for (int i = 0; i < AA.size(); i++)
        AA[i]=maxA;
    int FQ=0;
    while (true)
    {
        if (sz==-1) break;
        VecL [dl-1]++;
        if (VecL [dl-1] > VecM [dl-1])
        {
            ll = dl-1;
            while (ll>0)
            {
                VecL [ll] = 0;
                ll--;
                VecL [ll] ++;
                if (VecL[ll] <= VecM [ll]) break;
            }
        }
        Y.clear();
        f1=0;
        t2=0;
        for (int y=0; y<RofBfromA.size(); y++)
        {
            if (FQ==1) break;
            if (FQ==2) break;
            sx = 3*VecL[y];
            if (NumA[RofA[RofBfromA[y][sx]][RofBfromA[y][sx+1]]]==-1)
            {
                t2++;
                NumA[RofA[RofBfromA[y][sx]][RofBfromA[y][sx+1]]] = t2-1;
            }
            int x=0;
            AA[f1]=(RofA[RofBfromA[y][sx]][RofBfromA[y][sx+1]+x]);
            AA[f1+1]=(RofA[RofBfromA[y][sx]][RofBfromA[y][sx+1]+x+1]);
            f1+=2;
            if (npr)
                    Y.insert(std::pair <int, int> (AA[f1-2], AA[f1-1]));
            if (!npr)
                    Y.insert(std::pair <int, int> (std::min(AA[f1-2], AA[f1-1]), std::max(AA[f1-2], AA[f1-1])));
            if (Y.size()!=f1/2)
                {
                    if (VecL==VecM)
                    {
                        FQ = 2;
                        break;
                    }
                    for (int x1=RofBfromA.size()-1; x1>y; x1--)
                        VecL[x1]=VecM[x1];
                    if (VecL==VecM)
                    {
                        FQ = 2;
                        break;
                    }
                    for (int x1=0; x1<AA.size(); x1++)
                        NumA[AA[x1]]=-1;
                    FQ = 1;
                    break;
                }
                if (NumA[AA[f1-1]]==-1)
                {
                    t2++;
                    NumA[AA[f1-1]] = t2-1;
                }
                if (NumA[AA[f1-2]]>=MatrB.size())
                {
                    if (VecL==VecM)
                    {
                        FQ = 2;
                        break;
                    }
                    for (int x1=RofBfromA.size()-1; x1>y; x1--)
                        VecL[x1]=VecM[x1];
                    if (VecL==VecM)
                    {
                        FQ = 2;
                        break;
                    }
                    for (int x1=0; x1<AA.size(); x1++)
                        NumA[AA[x1]]=-1;
                    FQ = 1;
                    break;
                }
                if (NumA[AA[f1-1]]>=MatrB.size())
                {
                    if (VecL==VecM)
                    {
                        FQ = 2;
                        break;
                    }
                    for (int x1=RofBfromA.size()-1; x1>y; x1--)
                        VecL[x1]=VecM[x1];
                    if (VecL==VecM)
                    {
                        FQ = 2;
                        break;
                    }
                    for (int x1=0; x1<AA.size(); x1++)
                        NumA[AA[x1]]=-1;
                    FQ = 1;
                    break;
                }
                if (MatrB  [NumA[AA[f1-2]]][NumA[AA[f1-1]]]!=1)
                {
                    if (VecL==VecM)
                    {
                        FQ = 2;
                        break;
                    }
                    for (int x1=RofBfromA.size()-1; x1>y; x1--)
                        VecL[x1]=VecM[x1];
                    if (VecL==VecM)
                    {
                        FQ = 2;
                        break;
                    }
                    for (int x1=0; x1<AA.size(); x1++)
                        NumA[AA[x1]]=-1;
                    FQ = 1;
                    break;
                }
             }
        if (FQ==1)
        {
            FQ=0;
            continue;
        }
        if (FQ==2) break;
        if (BB.size()!=f1)
        {
            if (VecL==VecM)
            {
                break;
            }
            for (int x1=0; x1<AA.size(); x1++)
                NumA[AA[x1]]=-1;
            continue;
        }
        if (!npr)
        {
            for (int i=0; i<AA.size(); i+=2)
            {
                if (AA[i]>AA[i+1]) ZamVec(AA, i, i+1);
            }
        }
        RRR.clear();
        if (npr)
            for (int q=0; q<BB.size(); q+=2)
                for (int a=0; a<KratB[std::pair <int, int>(BB[q], BB[q+1]) ]; a++)
                    RRR.insert(std::pair<int, int> (AA[q],AA[q+1]));
        if (!npr)
            for (int q=0; q<BB.size(); q+=2)
                for (int a=0; a<KratB[std::pair <int, int>(std::min(BB[q], BB[q+1]), std::max(BB[q], BB[q+1])) ]; a++)
                    RRR.insert(std::pair<int, int> (AA[q],AA[q+1]));

        f3=0;
        for (auto it3=RRR.begin(); it3!=RRR.end(); it3++)
        {
            GR[f3]=(*it3).first;
            GR[f3+1]=(*it3).second;
            f3+=2;
        }
        if (sB!=f3)
        {
            if (VecL==VecM)
            {
                break;
            }
            for (int x1=0; x1<AA.size(); x1++)
                NumA[AA[x1]]=-1;
            continue;
        }
        if (GR.size()!=sB)
        {
            if (VecL==VecM)
            {
                break;
            }
            for (int x1=0; x1<AA.size(); x1++)
                NumA[AA[x1]]=-1;
            continue;
        }
        A_B.insert(GR);
        if (VecL==VecM)
        {
            break;
        }
        for (int x1=0; x1<AA.size(); x1++)
            NumA[AA[x1]]=-1;
    }
    double Time_ = time(NULL);
    TimeI = Time_ - Time;
    return A_B.size();
}


int main()
{
std::ifstream Chtenie;
std::ifstream ChtenieMetok;
std::vector <int> A;
std::vector <int> B;
std::map<int, std::string> Af_s = {};
std::map<int, std::string> Bf_s = {};

Chtenie.open("Graph_A.txt");            // чтение из файла
FileToVec (Chtenie, A);
Chtenie.close();
 Chtenie.open("Graph_B.txt");
FileToVec (Chtenie, B);
Chtenie.close();

ChtenieMetok.open ("Label_A.txt");      // чтение меток из файла
FileToMet (ChtenieMetok, Af_s);
ChtenieMetok.close();
ChtenieMetok.open ("Label_B.txt");
FileToMet (ChtenieMetok, Bf_s);
ChtenieMetok.close();

//    std::vector <int> A = {0, 1, 1, 2, 1, 3, 3, 1, 3, 4, 3, 5, 5, 3, 6, 7};  //пример 1
 //   std::vector <int> B = {0, 1, 0, 2, 2, 0, 3, 4};

    //std::vector <int> A = {0, 0, 0, 1, 1, 2, 1, 3,};                      //пример 2
    //std::vector <int> B = {0, 0, 0, 1};

 //std::vector <int> A = {0, 1, 1, 2, 1, 3, 3, 1, 3, 4, 3, 5, 5, 3};        //пример 3
 //std::vector <int> B = {0, 1, 0, 2, 2, 0};

    bool npr = true;
    int a = 0;

    std::set<std::vector <int>> A_B = * (new std::set<std::vector<int>>);

  //std::map<int, std::string> Af_s = {{1, "А"}};       //пример 1       метки А
   //std::map<int, std::string> Af_s = {{3, "А"}};     //пример 2
  //  std::map<int, std::string> Af_s = {};             //пример 3

    //std::map<int, std::string> Bf_s = {{0, "А"}};       // для примера 1 и 2   метки В
    //std::map<int, std::string> Bf_s = {};               // для примера 3


    a = PoiskGraph(A, B, A_B, npr, Af_s, Bf_s);        // поиск изоморфных

   std::cout<<a<<" - number of the finded subgraphs"<<std::endl;
    for (auto j=A_B.begin(); j!=A_B.end(); j++){
        for (int i=0; i<(*j).size(); i++){
    std::cout<<(*j)[i]<<"   ";}
    std::cout<<std::endl;}

std::ofstream Zapis;                                // запись в файл
Zapis.open("Result.txt", std::ios_base::out);
Zapis.clear();
Zapis.close();
Zapis.open("Result.txt", std::ios_base::app);
Zapis.clear();
for (auto j=A_B.begin(); j!=A_B.end(); j++)
    {
    VecToFile ((*j), Zapis);
    Zapis<<std::endl;
    }
Zapis.close();

std::cout<<std::endl<<TimeI;
}

