new;											@ 3 levels @		
t00=hsec/100;
load x0[]=start3.txt;							@ !!! starting design @
x0=reshape(x0,30,3);
L=maxc(x0[.,1]);								@ number of levels @
l1=L-1;
csn= rows(x0)/2 ;								@ number of choice sets @
b0=(ones(cols(x0),1)).*.(-1|0);
sgm= 1 ;
cn=cols(x0);
seed= 15 ;
s=2*cn;
ns=64;
load v[]=rndnorm.txt;
v=reshape(v,ns,s)';

nd= 2 ;											@ number of designs @
dd={};
ii=1;
do while ii<=nd;									@ start of loop @
t0=hsec/100;
xx=x0;											@ generates a random relabeling of x0 @
min=fcn(xx);
format /rd 6,3;
min;

/*  SWAPPING  */

m= cols(xx)-1 ;
ba= 2 ;
a={};
for i (1,ba^m-1,1);
	a=a|base(i)';
endfor;

label:
i=1;
do while i<=csn;
	q=1;
	do while q<=rows(a);
		x=xx;
		x[2*i-1:2*i,1:cols(xx)-1]=swv(x[2*i-1:2*i,1:cols(xx)-1],a[q,.]);
		fx=fcn(x);
		if fx<min;
			xx=x;
			min=fx;
			format /rd 3,0;"i=" i;
			format /rd 6,3;min;
			goto label;
		endif;
		q=q+1;
	endo;
	i=i+1;
endo;

/*  CYCLING  */

min=fcn(xx);

label2:
j=1;
do while j<=cols(xx);
	i=1;
	do while i<=csn;
		x=xx;
		x[2*i-1 2*i,j]=cy(x[2*i-1,j])|cy(x[2*i,j]);
		fx=fcn(x);
		if min>fx;
			min=fx;xx=x;
			format /rd 6,3;
			min;
			format /rd 3,0;"i=" i;"j=" j;
			goto label2;
		endif;
		x[2*i-1 2*i,j]=cy(x[2*i-1,j])|cy(x[2*i,j]);
		fx=fcn(x);
		if min>fx;
			min=fx;xx=x;
			format /rd 6,3;
			min;
			format /rd 3,0;"i=" i;"j=" j;
			goto label2;
		endif;
		x[2*i-1 2*i,j]=sw(xx[2*i-1 2*i,j],1);
		fx=fcn(x);
		if min>fx;
			min=fx;xx=x;
			format /rd 6,3;
			min;
			format /rd 3,0;"i=" i;"j=" j;
			goto label2;
		endif;
		x[2*i-1 2*i,j]=cy(x[2*i-1,j])|cy(x[2*i,j]);
		fx=fcn(x);
		if min>fx;
			min=fx;xx=x;
			format /rd 6,3;
			min;
			format /rd 3,0;"i=" i;"j=" j;
			goto label2;
		endif;
		x[2*i-1 2*i,j]=cy(x[2*i-1,j])|cy(x[2*i,j]);
		fx=fcn(x);
		if min>fx;
			min=fx;xx=x;
			format /rd 6,3;
			min;
			format /rd 3,0;"i=" i;"j=" j;
			goto label2;
		endif;
		i=i+1;
	endo;
	j=j+1;
endo;

dd=dd|xx;
"secs=" hsec/100-t0;
ii=ii+1;
endo;												@ end of loop @

@save baygre12=dd;@

"secs=" hsec/100-t00;

proc fcn(xx);
local b,edp,dx;
dx=dd|xx;
edp=0;
for i (1,ns,1);
	b=b0+sgm*v[.,i];
	edp=edp+dp(dx,b);
endfor;
retp(edp/ns);
endp;

proc dp(xx,b);
local mx,mm,mm1,s,i,t,rx,cx,d;
rx=rows(xx);
cx=cols(xx);
mx=reshape(xx,rx*cx,1);
mm=zeros(rx*cx,2);
i=1;
do while i<=rx*cx;
	if mx[i]==1;mm[i,1]=1;mm[i,2]=0;
	elseif mx[i]==2;mm[i,1]=0;mm[i,2]=1;
	else;mm[i,1]=-1;mm[i,2]=-1;
	endif;
	i=i+1;
endo;
mm1=reshape(mm,rx,l1*cx);
i=1;s=zeros(l1*cx,l1*cx);	
do while i<=rx/2;
	s=s+im(mm1[2*i-1 2*i,.],b);
	i=i+1;
endo;
d=det(s);
if d<=0;
	retp(exp(100));
else;
	retp(d^(-1/cols(s)));
endif;
endp;

proc im(z,b);
local nr,p;
nr=rows(z);
p=exp(z*b)/sumc(exp(z*b));
retp(z'*(diagrv(zeros(nr,nr),p)-p*p')*z);
endp;

proc sw(v,a);
local u;
if a==1;
	u=v[2,1];
	v[2,1]=v[1,1];
	v[1,1]=u;
endif;
retp(v);
endp;

proc swv(v,av);
local c,sv;
c=cols(v);
sv={};
for i (1,c,1);
	sv=sv~sw(v[.,i],av[i]);
endfor;
retp(sv);
endp;

proc cy(z);
if z==1;
z=2;
elseif z==2;
z=3;
else;
z=1;
endif;
retp(z);
endp;

proc perm(a,j);
if j==1;
	a=1*(a.<2)+2*(a.>2)+3*(a.>1).*(a.<3);
elseif j==2;
	a=1*(a.>2)+2*(a.>1).*(a.<3)+3*(a.<2);
elseif j==3;
	a=1*(a.>1).*(a.<3)+2*(a.<2)+3*(a.>2);
elseif j==4;
	a=1*(a.>1).*(a.<3)+2*(a.>2)+3*(a.<2);
elseif j==5;
	a=1*(a.>2)+2*(a.<2)+3*(a.>1).*(a.<3);
endif;
retp(a);
endp;

proc base(x);		@ gives a1,..,am as a column vector s.t. x=a1*b^(m-1)+...+a(m-1)*b+am @
local a,c,j;
a=zeros(m,1);
c=zeros(m,1);
a[m]=x%ba;
c[m]=a[m];
for j (m-1,1,-1);
	a[j]=((x-c[j+1])/ba^(m-j))%ba;
	c[j]=c[j+1]+a[j]*(ba^(m-j));
endfor;
retp(a);
endp;
