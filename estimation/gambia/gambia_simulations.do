
* Define stata output *

gen mean = 0
gen var = 0


* Start mata routine *

mata:

	T=8

	/* controls and other variables */

	ID = st_data(.,"id")

	y = st_data(., "y")

	z = (st_data(.,"year2"),st_data(.,"year3"),st_data(.,"year4"), st_data(.,"older"))


	x = (st_data(.,"treatment"),J(rows(z),1,1))

	CH = (st_data(.,"t"))

	TR = (st_data(.,"ever_treat"))



	/* First step: estimate parameters common across schools */

	deltaden=J(cols(z),cols(z),0)
	deltanum=J(cols(z),1,0)
	hz = J(cols(x),cols(z),0)

	j1=1
	maxID=max(ID)
	while (j1<=maxID){
		 y2 = select(y, (ID:==j1))
		 z2 = select(z, (ID:==j1))
		 x2 = (select(x, (ID:==j1)))
		 q2 = I(T)-x2*invsym(x2'x2)*x2' 
		 h2 = pinv(x2)

		 deltaden = deltaden+z2'q2*z2
		 deltanum = deltanum+z2'q2*y2
		 hz = hz+h2*z2
		 
		 j1=j1+1
		 }

	invdeltaden = invsym(deltaden)
	delta = invdeltaden*deltanum
	delta

	/* Second step: Estimate School Specific Effects Covariates Interacted with Fixed Effects */

	vardeltanum = J(cols(z),cols(z),0)
	meangamma = J(cols(x),1,0)
	varmeangamma1 = J(cols(x),cols(x),0)


	j1=1
	maxID=max(ID)
	while (j1<=maxID){
		 y2 = select(y, (ID:==j1))
		 z2 = select(z, (ID:==j1))
		 x2 = (select(x, (ID:==j1)))
		 q2 = I(T)-x2*invsym(x2'x2)*x2' 
		 h2 = pinv(x2)
		 
		 /* Standard errors common coefficients */
		 vardeltanum = vardeltanum+z2'*q2*(y2-z2*delta)*(y2-z2*delta)'*q2'*z2
		 
		 /* Overall mean individual effects */
		 meangamma = meangamma+h2*(y2-z2*delta)
		 
		 j1=j1+1
		}

	meangamma=meangamma/maxID
	MG=0
	st_view(MG,.,"mean")
	MG[.] = J(rows(y),1,meangamma[1,1])

	j1=1
	maxID=max(ID)
	while (j1<=maxID){
		 y2 = select(y, (ID:==j1))
		 z2 = select(z, (ID:==j1))
		 x2 = (select(x, (ID:==j1)))
		 q2 = I(T)-x2*invsym(x2'x2)*x2' 
		 h2 = pinv(x2)
		 
		 /* Standard errors average individual effects */
		 vect = h2*(y2-z2*delta-x2*meangamma)-hz*invdeltaden*z2'*q2*(y2-z2*delta)
		 varmeangamma1 = varmeangamma1+vect*vect'
			 
		 j1=j1+1
		}


	/* Standard Errors for Common coefficients */

	stddelta = sqrt(diagonal(invdeltaden*vardeltanum*invdeltaden'))
	display(" Coefficient -  Std. Dev. ")
	(delta,stddelta)

	/* Standard Error for Average of individual effects */

	stdmeangamma=sqrt(diagonal(varmeangamma1/maxID^2))

	varmeangamma1
	display(" Coefficient -  Std. Dev. ")
	(meangamma,stdmeangamma)

	/* Variances of Effects - Under Different Covariance Assumptions */

	/* Create Matrices */
	invxx=J(cols(x),cols(x),0)
	sigma2u=0
	vargamma=J(cols(x),cols(x),0)
	varsigma2u=0
	varvargammaiid=J(2*cols(x)^2+1,2*cols(x)^2+1,0)
	invxxrobust=J(cols(x),cols(x),0)
	varvargammarobust=J(cols(x)^2,cols(x)^2,0)
	varutvL = J(T^2,1,0)
	varvarutvL=J(T^2,T^2,0)
	varutvgL = J(cols(x)^2,1,0)
	varvargammatvL=J(cols(x)^2,cols(x)^2,0)
	varbeta=J(cols(x)-1,cols(x)-1,0)
	varvarbeta=J((cols(x)-1)^2,(cols(x)-1)^2,0)


	j1=1
	maxID=max(ID)
	while (j1<=maxID){
		 y2 = select(y, (ID:==j1))
		 z2 = select(z, (ID:==j1))
		 x2 = (select(x, (ID:==j1)))
		 q2 = I(T)-x2*invsym(x2'x2)*x2' 
		 h2 = pinv(x2)
		 q2b = I(T)-J(T,T,1)/T
		 h2b=pinv(q2b*x2[.,1])
		 
		 /* Variance errors, iid */
		 sigma2u = sigma2u+(y2-z2*delta)'q2*(y2-z2*delta)/(T-cols(x))

		 /* Variance individual effects, iid */
		 vargamma = vargamma+h2*(y2-z2*delta-x2*meangamma)*(y2-z2*delta-x2*meangamma)'*h2'
		 invxx = invxx+invsym(x2'x2)

		 /* Variance individual effects, robust */
		 invxxrobust = invxxrobust+(y2-z2*delta)'q2*(y2-z2*delta)/(T-cols(x))*invsym(x2'x2)

		 /* Variance beta in first differences, robust */
		 varbeta = varbeta+h2b*q2b*(y2-z2*delta-x2*meangamma)*(y2-z2*delta-x2*meangamma)'*q2b'*h2b'/*
		 */ -(y2-z2*delta)'q2*(y2-z2*delta)/(T-cols(x))*h2b*q2b*h2b'

		 j1=j1+1
		}


	j1=1
	maxID=max(ID)
	while (j1<=maxID){
		 y2 = select(y, (ID:==j1))
		 z2 = select(z, (ID:==j1))
		 x2 = (select(x, (ID:==j1)))
		 q2 = I(T)-x2*invsym(x2'x2)*x2' 
		 h2 = pinv(x2)
		 q2b = I(T)-J(T,T,1)/T
		 h2b=pinv(q2b*x2[.,1])


		 /* std variance errors, iid */
		 varsigma2u = varsigma2u+((y2-z2*delta)'q2*(y2-z2*delta)/(T-cols(x)))^2

		 /* std variance individual effects, iid */
		 vect=vec(h2*(y2-z2*delta-x2*meangamma)*(y2-z2*delta-x2*meangamma)'*h2')
		 vect=vect\(y2-z2*delta)'q2*(y2-z2*delta)/(T-cols(x))
		 vect=vect\vec(invsym(x2'x2))
		 varvargammaiid = varvargammaiid+vect*vect'

		 /* std variance individual effects, robust */
		 vect=vec(h2*(y2-z2*delta-x2*meangamma)*(y2-z2*delta-x2*meangamma)'*h2'/*
				*/-(y2-z2*delta)'q2*(y2-z2*delta)/(T-cols(x))*invsym(x2'x2))
		 varvargammarobust = varvargammarobust+vect*vect'

		 /* std variance beta in first differences, robust */
		 vect = h2b*q2b*(y2-z2*delta-x2*meangamma)*(y2-z2*delta-x2*meangamma)'*q2b'*h2b'/*
		 */ -(y2-z2*delta)'q2*(y2-z2*delta)/(T-cols(x))*h2b*q2b*h2b'
		   varvarbeta = varvarbeta+vect*vect'
		 

		 j1=j1+1
		}




	/* Variance errors, iid errors */
	sigma2u=sigma2u/maxID
	stdsigma2u=sqrt(varsigma2u/maxID^2)
	display(" Coefficient -  Std. Dev. ")
	sigma2u,stdsigma2u

	/* Variance individual effects, iid errors */
	vargamma=vargamma/maxID
	vargammaiid=vargamma-sigma2u*invxx/maxID
	M=I(cols(x)^2),(-vec(invxx/maxID)),(-sigma2u*I(cols(x)^2))
	stdvargammaiid = sqrt(diagonal(M*varvargammaiid*M'/maxID^2-vec(vargammaiid)*vec(vargammaiid)'/maxID))
	display(" Coefficient -  Std. Dev. ")
	vec(vargammaiid),stdvargammaiid

	/* Variance individual effects, robust */
	vargammarobust = vargamma-invxxrobust/maxID
	stdvargammarobust = sqrt(diagonal(varvargammarobust/maxID^2-vec(vargammarobust)*vec(vargammarobust)'/maxID))
	display(" Coefficient -  Std. Dev. ")
	vec(vargammarobust),stdvargammarobust

	/* Export Prefered Variance To Stata (out of mata) */
	V=0
	st_view(V,.,"var")
	V[.] = J(rows(y),1,vargammarobust[1,1])
end
