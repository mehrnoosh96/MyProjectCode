sets
i nodes /1*5/
alias (i,j,k,m,p,s) ;

parameter
alfa discount factor for line-haul movement between hubs /0.2/
pp the number of hubs to locate /2/

teta(i)
/
1  200
2  149
3  40
4  89
5  190/ ;

table h(i,k) demand between origin i to destination j
        1     2    3    4    5
1       0    32   25   77   20

2       22    0   36   28   42

3       36   36   0   30    65

4       67    28   30   0   42

5       15    34   56   46   0   ;


table c(i,j) unit cost of local (non-hub to hub) movement between nodes i to j
        1       2      3    4     5

  1     0       27    23    19    15

  2     34      0     26    24    32

  3     21      18    0     25    49

  4     27      24     25    0    26

  5     14      29     53    38    0;



variable
zz total cost;
positive variable
z(i,j,k,m);

binary variables
x(j) a hub is located at node j
y(i,j) node i is connected to a hub located at node j;



equation
obj
co1
co2
co4
co5
*co6
co7


;

obj..           zz=e=sum((i,j,k,m),(c(i,k)+c(k,m)+c(m,j))*h(i,j)*z(i,j,k,m));
co1..                                 sum(k,x(k))=e=pp   ;
co2(i,j)$(ord(i)<>ord(j))..           sum((k,m),z(i,j,k,m))=e=1;
co4(i,j,k,m)$(ord(i)<>ord(j)) ..      z(i,j,k,m)=l=x(m);
co5(i,j,k,m)$(ord(i)<>ord(j)) ..      z(i,j,k,m)=l=x(k);
*co6(i,j,k,m) $(ord(i)<>ord(j))..      y(i,k)+y(j,m)-2*z(i,j,k,m)=g=0;
co7(i,k).. sum((m,j),h(i,j)*z(i,j,k,m))+sum((p,s),h(p,i)*z(p,i,s,k))=g=teta(k)*y(i,k);





model PHL /all/;
option mip=cplex;
option optca=0;
option optcr=0;
solve PHL using mip min zz;
display x.l,z.l,zz.l;

