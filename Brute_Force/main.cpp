#include <fstream>
using std::ifstream;
#include <algorithm>
using std::is_sorted;
#include <assert.h>
#include <iterator>
using std::reverse_iterator;
#include <iostream>
using std::ostream;
using std::cout;
using std::endl;
#include <numeric>
using std::iota;
#include <set>
using std::set;
#include <utility>
using std::begin;
using std::end;
#include <vector>
using std::vector;
using std::pair;
using std::string;

int colv = 3, colr = 2, colv_G_A = 8, colr_G_A = 5;      // кол-во вершин и ребер ддя графов В и А  начиная с 0
vector<vector<int>> Na_1;
vector<vector<int>> New_p;

I()
{
    Na_1.resize(colv_G_A);
    for (int i=0; i<colv_G_A; i++)
    {
        (Na_1[i]).resize(colv_G_A);
    }
    New_p.resize(colv_G_A);
    for (int i=0; i<colv_G_A; i++)
    {
        (New_p[i]).resize(colv_G_A);
    }
}

int Matrica (string Name_1, int n, int m)
{
    int s, p;
    ifstream in(Name_1);

    vector<vector<int>> Na;
    Na.resize(n);
    for (int i=0; i<n; i++)
    {
        (Na[i]).resize(n);
    }
    for (int k = 0; k<m; k++)
    {
        in >> s >> p;
        Na [s] [p] = 1;                 //нумерация с 1 - то -1 ДОБАВИТЬ к s и p, ИНАЧЕ УБРАТЬ.
    }

    for (int i = 0; i<n; i++)
    {
        for (int j = 0; j<n; j++)
        {
            Na_1 [i] [j] = Na [i] [j];
        }
    }
    in.close();
    return 0;
};

int ISO(vector<vector<int>> AB_L,  vector<vector<int>> AB_P, int RM)
{
    vector <int> p;
    p.resize(RM);
    for(int k=0; k<RM; k++)
    {
        p[k] = k;
    };
    int x,AA = 0 ;
    do
    {
        x=0;
        for (int i=0; i<RM; i++)
        {
            for (int j=0; j<RM; j++)
            {
                if ( AB_L [i] [j] ==  AB_P [ p [i] ][ p [j] ])
                {
                    x++;
                }
            }
        }
        if (x==(RM*RM))
        {
            cout<<"izo"<<x<<endl;
        }
    }
    while (std::next_permutation(p.begin(), p.end()));
    return 0;
}

int Matrica_ (vector <int> Vec, int n)
{
    int VS =  Vec.size();
    int s = 0,p = 0;
    vector<vector<int>> Na;
    Na.resize(n);
    for (int i=0; i<n; i++)
    {
        (Na[i]).resize(n);
    }
    for (int k = 0; k<VS; k+=2)
    {
        s = Vec[k];
        p = Vec[k+1];
        Na [s] [p] = 1;                     //нумерация с 1 - то -1 ДОБАВИТЬ к s и p, ИНАЧЕ УБРАТЬ.
    }
    for (int i = 0; i<n; i++)
    {
        for (int j = 0; j<n; j++)
        {
            New_p [i][j] =  Na[i][j];
        }
    }
    return 0;
};

template< class It >
auto ascend_ordered( const int n_digits, const It begin, const It end )
-> bool
{
    using R_it = reverse_iterator< It >;
    const R_it r_begin  = R_it( end );
    const R_it r_end    = R_it( begin );

    int max_digit = n_digits - 1;
    for( R_it it = r_begin ; it != r_end; ++it )
    {
        if( *it < max_digit )
        {
            ++*it;
            const int n_further_items = it - r_begin;
            for( It it2 = end - n_further_items; it2 != end; ++it2 )
            {
                *it2 = *(it2 - 1) + 1;
            }
            return true;
        }
        --max_digit;
    }
    return false;
}

int main()
{
    I();
    Matrica ("Graph_B.txt", colv_G_A, colr);
    ifstream in("Graph_A.txt");
    vector<pair<int,int>> V_;
    int ab;
    V_.resize(colr_G_A);
    for (int i = 0; i<colr_G_A*2; i+=2)
    {
        in>>V_[i/2].first;
        in>>V_[i/2].second;
    }
    sort(V_.begin(), V_.end());
    assert( is_sorted( begin( V_ ), end( V_ ) ) );
    const int k = colr;
    const int n = V_.size();
    vector<int> indices( k );
    iota( indices.begin(), indices.end(), 0 );
    set<vector<int>> encountered;
    int ii_ =0;
    vector<int> ac;
    for( ;; )
    {
        vector<pair<int,int>> current;
        for( int const i : indices )
        {
            current.push_back(V_[i] );
        }
        {
            for (int i=0; i<current.size(); i++)
            {
                ac.push_back(current[i].first);
                ac.push_back(current[i].second);
            }
            for (int i=0; i<ac.size(); i++)
            {
            }
            encountered.insert( ac );
            ac.clear();
        }
        if( not ascend_ordered( n, begin( indices ), end( indices ) ) )
        {
            break;
        }
    }

    for (auto V =  encountered.begin(); V!=encountered.end(); V++)
    {

        Matrica_ ((*V), /*мб*/colv_G_A);
        ISO(New_p, Na_1, colv_G_A);
    }
    return 0;
}
