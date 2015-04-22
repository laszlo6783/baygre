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
transf={3 3 3};

IM0=zeros(s*ns,s*csn);	/*	Az aktualis choice-set-ek IM-je	*/
IMTot=zeros(s*ns,s);		/*	Az elozo designok osszesitett IM-je	*/
IMTotA=zeros(s*ns,s);		/*	Az aktualis osszesitett IM	*/
IMTemp=zeros(s*ns,s);		/*	A temporalis osszesitett IM	*/
IM0Temp=zeros(s*ns,s);	/*	A temporalis choice-set IM-je	*/
transf={3 3 3};

for i (1,ns,1);
for ii (1,csn,1);
	b=b0+sgm*v[.,i];
	IM0[6*i-5:6*i,6*ii-5:6*ii]=improc(xtoz(X0[2*ii-1:2*ii,.],transf),b);
endfor;
endfor;

for i (1,ns,1);
for ii (1,csn,1);
	IMTot[6*i-5:6*i,.]=IMTot[6*i-5:6*i,.]+IM0[6*i-5:6*i,6*ii-5:6*ii];
endfor;
endfor;

IMT00=IMTot;	/*legelso, ha UGYANAZT a designt hasznaljuk mindig, s nem keverjuk meg, akkor hasznalhato */

nd= 4 ;											@ number of designs @
dd={};
ii=1;
do while ii<=nd;									@ start of loop @
t0=hsec/100;
xx=X0;											@ generates a random relabeling of x0 @
min=fcnN(IMTot);
format /rd 6,3;
min;

if ii==1;
IMTot=zeros(s*ns,s);
endif;
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
		{fx,IMTemp}=fcn2(IMTot,x,v);
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

{min,IMTemp}=fcn2(IMTot,xx,v);

label2:
j=1;
do while j<=cols(xx);
	i=1;
	do while i<=csn;
		x=xx;
		x[2*i-1 2*i,j]=cy(x[2*i-1,j])|cy(x[2*i,j]);
		{fx,IMTemp}=fcn2(IMTot,x,v);
		if min>fx;
			min=fx;xx=x;
			format /rd 6,3;
			min;
			format /rd 3,0;"i=" i;"j=" j;
			goto label2;
		endif;
		x[2*i-1 2*i,j]=cy(x[2*i-1,j])|cy(x[2*i,j]);
		{fx,IMTemp}=fcn2(IMTot,x,v);
		if min>fx;
			min=fx;xx=x;
			format /rd 6,3;
			min;
			format /rd 3,0;"i=" i;"j=" j;
			goto label2;
		endif;
		x[2*i-1 2*i,j]=sw(xx[2*i-1 2*i,j],1);
		{fx,IMTemp}=fcn2(IMTot,x,v);
		if min>fx;
			min=fx;xx=x;
			format /rd 6,3;
			min;
			format /rd 3,0;"i=" i;"j=" j;
			goto label2;
		endif;
		x[2*i-1 2*i,j]=cy(x[2*i-1,j])|cy(x[2*i,j]);
		{fx,IMTemp}=fcn2(IMTot,x,v);
		if min>fx;
			min=fx;xx=x;
			format /rd 6,3;
			min;
			format /rd 3,0;"i=" i;"j=" j;
			goto label2;
		endif;
		x[2*i-1 2*i,j]=cy(x[2*i-1,j])|cy(x[2*i,j]);
		{fx,IMTemp}=fcn2(IMTot,x,v);
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
IMTot=IMTemp;
dd=dd|xx;
"secs=" hsec/100-t0;
ii=ii+1;
endo;												@ end of loop @

@save baygre12=dd;@

"secs=" hsec/100-t00;

PROC(3) = fcn3(PIM,x,LL,v);
LOCAL ;
ENDP;

PROC(2) = fcn2(PIM,actx,v);
LOCAL j,jj,derr,IMTemp;
IM0=zeros(s*ns,s*csn);
IMTemp=zeros(s*ns,s);
for j (1,ns,1);
for jj (1,csn,1);
	b=b0+sgm*v[.,j];
	IM0[6*j-5:6*j,6*jj-5:6*jj]=improc(xtoz(actx[2*jj-1:2*jj,.],transf),b);
endfor;
endfor;

for j (1,ns,1);
for jj (1,csn,1);
	IMTemp[6*j-5:6*j,.]=IMTemp[6*j-5:6*j,.]+IM0[6*j-5:6*j,6*jj-5:6*jj];
endfor;
endfor;
IMTemp=PIM+IMTemp;
derr=fcnN(IMtemp);
RETP(derr,IMTemp);
ENDP;

PROC(1)	= fcnN(imt);
LOCAL j,d,derr;
derr=0;
for j (1,ns,1);
	d=det(imt[6*j-5:6*j,.]);
	if d<0;
		derr=derr+exp(100);
	else;
		derr=derr+d^(-1/6);
	endif;
endfor;
RETP(derr/ns);
ENDP;

PROC(1)  = improc(Z,b);
	LOCAL col,row,j,imm;
	col=cols(Z);
	row = rows(Z);
	j = 1;
	imm = zeros(col,col);
	DO WHILE j<=row/2;
	    imm=imm+im(Z[2*j-1 2*j,.],b);
	    j=j+1;
	ENDO;
	RETP(imm);
ENDP;

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

/*	This procedure transforms matrix X in coded matrix Z*/
/*	levels is a vector with length of parameter vector and value of how many levels we have*/
PROC(1) = xtoz (inX,levels);
	LOCAL c,r,c2,outZ,k;

	c=cols(inX);
	r=rows(inX);
	c2=cols(levels);

	if c ne c2;
		print "Size of the parameters differ";;
	endif;

	for i(1,c2,1);
		if levels[i]==3;
			c2=c2+1;
		endif;
	endfor;
	
	outZ=zeros(r,c2);
	
	k=1; /*	Column pozition in Z matrix	*/

	for i(1,c,1);
		if levels[i]==2;
			for row(1,r,1);
				if inX[row,i]==1;
					outZ(row,k)=-1;
				else;
					outZ(row,k)=1;
				endif;
			endfor;
			k=k+1;
		else;
			for row(1,r,1);
				if inX[row,i] == 1;
					outZ[row,k]=1;
					outZ[row,k+1]=0;
				elseif inX[row,i] == 2;
					outZ[row,k]=0;
					outZ[row,k+1]=1;
				else;
					outZ[row,k]=-1;
					outZ[row,k+1]=-1;
				endif;
			endfor;
			k=k+2;
		endif;
	endfor;

	RETP (outZ);

ENDP;
