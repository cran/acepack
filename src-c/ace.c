/* ace.f -- translated by f2c (version 19950110).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct prams_1_ {
    integer itape, maxit, nterm;
    doublereal span, alpha, big;
};

#define prams_1 (*(struct prams_1_ *) &prams_)

/* Initialized data */

struct {
    integer e_1[3];
    doublereal e_2[3];
    } prams_ = { -6, 20, 3, 0., 0., 1e20 };


/* Table of constant values */

static integer c__1 = 1;

/*     real -> double precision conversion for R use */
/*     <TSL> */
/*     mortran 2.0     (version of 6/24/75) */
/* Subroutine */ int mace_(p, n, x, y, w, l, delrsq, ns, tx, ty, rsq, ierr, m,
	 z)
integer *p, *n;
doublereal *x, *y, *w;
integer *l;
doublereal *delrsq;
integer *ns;
doublereal *tx, *ty, *rsq;
integer *ierr, *m;
doublereal *z;
{
    /* Format strings */
    static char fmt_670[] = "(\002 ierr=6: l(\002i2,\002) =\002g12.4,\002 mu\
st be in the range (-5, 5).\002)";
    static char fmt_650[] = "(\002 ierr=4: l(\002i2,\002) must be nonzero\
.\002)";
    static char fmt_660[] = "(\002 ierr=5: at least one l(1)-l(\002i2,\002) \
must be nonzero.\002)";
    static char fmt_620[] = "(\002 ierr=1: sum of weights (w) not positive\
.\002)";
    static char fmt_590[] = "(\0020eigensolution \002i2,\002:\002)";
    static char fmt_630[] = "(\002 ierr=2: y has zero variance.\002)";
    static char fmt_640[] = "(\002 ierr=3: ty(.,\002i2,\002) has zero varian\
ce.\002)";
    static char fmt_610[] = "(\002     iteration \002i2,\002   r**2  =  1 - \
e**2  =\002g12.4)";
    static char fmt_600[] = "(\002 eigensolution \002i2,\002   r**2  =  1 - \
e**2  =\002g12.4)";

    /* System generated locals */
    integer m_dim1, m_offset, x_dim1, x_offset, ty_dim1, ty_offset, tx_dim1, 
	    tx_dim2, tx_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();
    double sqrt();

    /* Local variables */
    static integer iter;
    static doublereal rsqi;
    extern /* Subroutine */ int sort_();
    static integer i, j, k;
    extern /* Subroutine */ int scale_();
    static doublereal ct[10];
    static integer is, js, np;
    static doublereal sm;
    static integer nt;
    static doublereal sv, sw;
    extern /* Subroutine */ int smothr_();
    static integer pp1;
    static doublereal sw1, cmn, cmx;
    static integer nit, ism1;

    /* Fortran I/O blocks */
    static cilist io___7 = { 0, 0, 0, fmt_670, 0 };
    static cilist io___8 = { 0, 0, 0, fmt_650, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_660, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_620, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_590, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_620, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_630, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_640, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_640, 0 };
    static cilist io___27 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_600, 0 };



/*   subroutine mace(p,n,x,y,w,l,delrsq,ns,tx,ty,rsq,ierr,m,z) */
/* ------------------------------------------------------------------ */

/* estimate multiple optimal transformations for regression and */
/* correlation by alternating conditional expectation estimates. */

/* version 3/28/85. */

/* breiman and friedman, journal of the american statistical */
/* association (september, 1985) */

/* coded  and copywrite (c) 1985 by: */

/*                        jerome h. friedman */
/*                     department of statistics */
/*                               and */
/*                stanford linear accelerator center */
/*                        stanford university */

/* all rights reserved. */


/* input: */

/*    n : number of observations. */
/*    p : number of predictor variables for each observation. */
/*    x(p,n) : predictor data matrix. */
/*    y(n) : response values for the observations. */
/*       missing values are signified by a value (response or */
/*       predictor) greater than or equal to big. */
/*       (see below - default, big = 1.0e20) */
/*    w(n) : weights for the observations. */
/*    l(p+1) : flag for each variable. */
/*       l(1) through l(p) : predictor variables. */
/*       l(p+1) : response variable. */
/*       l(i)=0 => ith variable not to be used. */
/*       l(i)=1 => ith variable assumes orderable values. */
/*       l(i)=2 => ith variable assumes circular (periodic) values */
/*                 in the range (0.0,1.0) with period 1.0. */
/*       l(i)=3 => ith variable transformation is to be monotone. */
/*       l(i)=4 => ith variable transformation is to be linear. */
/*       l(i)=5 => ith variable assumes categorical (unorderable) values. 
*/
/*   delrsq : termination threshold. iteration stops when */
/*       rsq changes less than delrsq in nterm */
/*       consecutive iterations (see below - default, nterm=3). */
/*   ns : number of eigensolutions (sets of transformations). */

/* output: */

/*   tx(n,p,ns) : predictor transformations. */
/*      tx(j,i,k) = transformed value of ith predictor for jth obs */
/*                  for kth eigensolution. */
/*   ty(n,ns) = response transformations. */
/*      ty(j,k) = transformed response value for jth observation */
/*                for kth eigensolution. */
/*   rsq(ns) = fraction of variance(ty<y>) */
/*                       p */
/*         explained by sum tx(i)<x(i)>  for each eigensolution. */
/*                      i=1 */
/*   ierr : error flag. */
/*      ierr = 0 : no errors detected. */
/*      ierr > 0 : error detected - see format statements below. */

/* scratch: */

/*    m(n,p+1), z(n,12) : internal working storage. */

/* note: mace uses an iterative procedure for solving the optimization */
/*    problem. default starting transformations are ty(j,k)=y(j), */
/*    tx(j,i,k)=x(i,j) : j=1,n, i=1,p, k=1,ns. other starting transformat 
*/
/*    can be specified (if desired) for either the response and/or any of 
*/
/*    the predictor variables. this is signaled by negating the */
/*    corresponding l(i) value and storing the starting transformed */
/*    values in the corresponding array (ty(j,k), tx(j,i,k)) before */
/*    calling mace. */

/* ------------------------------------------------------------------ */

    /* Parameter adjustments */
    z_dim1 = *n;
    z_offset = z_dim1 + 1;
    z -= z_offset;
    m_dim1 = *n;
    m_offset = m_dim1 + 1;
    m -= m_offset;
    --w;
    --y;
    x_dim1 = *p;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    --l;
    --rsq;
    ty_dim1 = *n;
    ty_offset = ty_dim1 + 1;
    ty -= ty_offset;
    tx_dim1 = *n;
    tx_dim2 = *p;
    tx_offset = tx_dim1 * (tx_dim2 + 1) + 1;
    tx -= tx_offset;

    /* Function Body */
    *ierr = 0;
    pp1 = *p + 1;
    sm = (float)0.;
    sv = sm;
    sw = sv;
    sw1 = sw;
    i__1 = pp1;
    for (i = 1; i <= i__1; ++i) {
	if (l[i] >= -5 && l[i] <= 5) {
	    goto L10;
	}
	*ierr = 6;
	if (prams_1.itape > 0) {
	    io___7.ciunit = prams_1.itape;
	    s_wsfe(&io___7);
	    do_fio(&c__1, (char *)&i, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&l[i], (ftnlen)sizeof(integer));
	    e_wsfe();
	}
L10:
	;
    }
    if (*ierr != 0) {
	return 0;
    }
    if (l[pp1] != 0) {
	goto L20;
    }
    *ierr = 4;
    if (prams_1.itape > 0) {
	io___8.ciunit = prams_1.itape;
	s_wsfe(&io___8);
	do_fio(&c__1, (char *)&pp1, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    return 0;
L20:
    np = 0;
    i__1 = *p;
    for (i = 1; i <= i__1; ++i) {
	if (l[i] != 0) {
	    ++np;
	}
/* L30: */
    }
    if (np > 0) {
	goto L40;
    }
    *ierr = 5;
    if (prams_1.itape > 0) {
	io___10.ciunit = prams_1.itape;
	s_wsfe(&io___10);
	do_fio(&c__1, (char *)&(*p), (ftnlen)sizeof(integer));
	e_wsfe();
    }
    return 0;
L40:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	sw += w[j];
/* L50: */
    }
    if (sw > (float)0.) {
	goto L60;
    }
    *ierr = 1;
    if (prams_1.itape > 0) {
	io___12.ciunit = prams_1.itape;
	s_wsfe(&io___12);
	e_wsfe();
    }
    return 0;
L60:
    i__1 = *ns;
    for (is = 1; is <= i__1; ++is) {
	if (prams_1.itape > 0) {
	    io___14.ciunit = prams_1.itape;
	    s_wsfe(&io___14);
	    do_fio(&c__1, (char *)&is, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    if (l[pp1] > 0) {
		ty[j + is * ty_dim1] = y[j];
	    }
/* L70: */
	}
	i__2 = *p;
	for (i = 1; i <= i__2; ++i) {
	    if (l[i] != 0) {
		goto L90;
	    }
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		tx[j + (i + is * tx_dim2) * tx_dim1] = (float)0.;
/* L80: */
	    }
	    goto L170;
L90:
	    if (l[i] <= 0) {
		goto L110;
	    }
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		tx[j + (i + is * tx_dim2) * tx_dim1] = x[i + j * x_dim1];
/* L100: */
	    }
L110:
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		if (tx[j + (i + is * tx_dim2) * tx_dim1] >= prams_1.big) {
		    goto L120;
		}
		sm += w[j] * tx[j + (i + is * tx_dim2) * tx_dim1];
		sw1 += w[j];
L120:
		;
	    }
	    if (sw1 > (float)0.) {
		goto L140;
	    }
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		tx[j + (i + is * tx_dim2) * tx_dim1] = (float)0.;
/* L130: */
	    }
	    sm = (float)0.;
	    sw1 = sm;
	    goto L170;
L140:
	    sm /= sw1;
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		if (tx[j + (i + is * tx_dim2) * tx_dim1] >= prams_1.big) {
		    goto L150;
		}
		tx[j + (i + is * tx_dim2) * tx_dim1] -= sm;
		goto L160;
L150:
		tx[j + (i + is * tx_dim2) * tx_dim1] = (float)0.;
L160:
		;
	    }
	    sm = (float)0.;
	    sw1 = sm;
L170:
	    ;
	}
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    if (ty[j + is * ty_dim1] >= prams_1.big) {
		goto L180;
	    }
	    sm += w[j] * ty[j + is * ty_dim1];
	    sw1 += w[j];
L180:
	    ;
	}
	if (sw1 > (float)0.) {
	    goto L190;
	}
	*ierr = 1;
	if (prams_1.itape > 0) {
	    io___15.ciunit = prams_1.itape;
	    s_wsfe(&io___15);
	    e_wsfe();
	}
	return 0;
L190:
	sm /= sw1;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    if (ty[j + is * ty_dim1] >= prams_1.big) {
		goto L200;
	    }
	    ty[j + is * ty_dim1] -= sm;
	    goto L210;
L200:
	    ty[j + is * ty_dim1] = (float)0.;
L210:
	    ;
	}
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
/* Computing 2nd power */
	    d__1 = ty[j + is * ty_dim1];
	    sv += w[j] * (d__1 * d__1);
/* L220: */
	}
	sv /= sw;
	if (sv <= (float)0.) {
	    goto L230;
	}
	sv = (float)1. / sqrt(sv);
	goto L260;
L230:
	if (l[pp1] <= 0) {
	    goto L240;
	}
	*ierr = 2;
	if (prams_1.itape > 0) {
	    io___16.ciunit = prams_1.itape;
	    s_wsfe(&io___16);
	    e_wsfe();
	}
	goto L250;
L240:
	*ierr = 3;
	if (prams_1.itape > 0) {
	    io___17.ciunit = prams_1.itape;
	    s_wsfe(&io___17);
	    do_fio(&c__1, (char *)&is, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
L250:
	return 0;
L260:
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    ty[j + is * ty_dim1] *= sv;
/* L270: */
	}
	if (is != 1) {
	    goto L310;
	}
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    m[j + pp1 * m_dim1] = j;
	    z[j + (z_dim1 << 1)] = y[j];
/* L280: */
	}
	sort_(&z[(z_dim1 << 1) + 1], &m[pp1 * m_dim1 + 1], &c__1, n);
	i__2 = *p;
	for (i = 1; i <= i__2; ++i) {
	    if (l[i] == 0) {
		goto L300;
	    }
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		m[j + i * m_dim1] = j;
		z[j + (z_dim1 << 1)] = x[i + j * x_dim1];
/* L290: */
	    }
	    sort_(&z[(z_dim1 << 1) + 1], &m[i * m_dim1 + 1], &c__1, n);
L300:
	    ;
	}
L310:
	scale_(p, n, &w[1], &sw, &ty[is * ty_dim1 + 1], &tx[(is * tx_dim2 + 1)
		 * tx_dim1 + 1], delrsq, p, &z[z_dim1 * 5 + 1], &z[z_dim1 * 6 
		+ 1]);
	rsq[is] = (float)0.;
	iter = 0;
	prams_1.nterm = min(prams_1.nterm,10);
	nt = 0;
	i__2 = prams_1.nterm;
	for (i = 1; i <= i__2; ++i) {
	    ct[i - 1] = (float)100.;
/* L320: */
	}
L330:
	++iter;
	nit = 0;
L340:
	rsqi = rsq[is];
	++nit;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    z[j + z_dim1 * 5] = ty[j + is * ty_dim1];
	    i__3 = *p;
	    for (i = 1; i <= i__3; ++i) {
		if (l[i] != 0) {
		    z[j + z_dim1 * 5] -= tx[j + (i + is * tx_dim2) * tx_dim1];
		}
/* L350: */
	    }
/* L360: */
	}
	i__2 = *p;
	for (i = 1; i <= i__2; ++i) {
	    if (l[i] == 0) {
		goto L420;
	    }
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		k = m[j + i * m_dim1];
		z[j + z_dim1] = z[k + z_dim1 * 5] + tx[k + (i + is * tx_dim2) 
			* tx_dim1];
		z[j + (z_dim1 << 1)] = x[i + k * x_dim1];
		z[j + (z_dim1 << 2)] = w[k];
/* L370: */
	    }
	    i__4 = (i__3 = l[i], abs(i__3));
	    smothr_(&i__4, n, &z[(z_dim1 << 1) + 1], &z[z_offset], &z[(z_dim1 
		    << 2) + 1], &z[z_dim1 * 3 + 1], &z[z_dim1 * 6 + 1]);
	    sm = (float)0.;
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		sm += z[j + (z_dim1 << 2)] * z[j + z_dim1 * 3];
/* L380: */
	    }
	    sm /= sw;
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		z[j + z_dim1 * 3] -= sm;
/* L390: */
	    }
	    sv = (float)0.;
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
/* Computing 2nd power */
		d__1 = z[j + z_dim1] - z[j + z_dim1 * 3];
		sv += z[j + (z_dim1 << 2)] * (d__1 * d__1);
/* L400: */
	    }
	    sv = (float)1. - sv / sw;
	    if (sv <= rsq[is]) {
		goto L420;
	    }
	    rsq[is] = sv;
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		k = m[j + i * m_dim1];
		tx[k + (i + is * tx_dim2) * tx_dim1] = z[j + z_dim1 * 3];
		z[k + z_dim1 * 5] = z[j + z_dim1] - z[j + z_dim1 * 3];
/* L410: */
	    }
L420:
	    ;
	}
	if (np == 1 || rsq[is] - rsqi <= *delrsq || nit >= prams_1.maxit) {
	    goto L430;
	}
	goto L340;
L430:
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    k = m[j + pp1 * m_dim1];
	    z[j + (z_dim1 << 1)] = y[k];
	    z[j + (z_dim1 << 2)] = w[k];
	    z[j + z_dim1] = (float)0.;
	    i__3 = *p;
	    for (i = 1; i <= i__3; ++i) {
		if (l[i] != 0) {
		    z[j + z_dim1] += tx[k + (i + is * tx_dim2) * tx_dim1];
		}
/* L440: */
	    }
/* L450: */
	}
	i__3 = (i__2 = l[pp1], abs(i__2));
	smothr_(&i__3, n, &z[(z_dim1 << 1) + 1], &z[z_offset], &z[(z_dim1 << 
		2) + 1], &z[z_dim1 * 3 + 1], &z[z_dim1 * 6 + 1]);
	if (is <= 1) {
	    goto L490;
	}
	ism1 = is - 1;
	i__2 = ism1;
	for (js = 1; js <= i__2; ++js) {
	    sm = (float)0.;
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		k = m[j + pp1 * m_dim1];
		sm += w[k] * z[j + z_dim1 * 3] * ty[k + js * ty_dim1];
/* L460: */
	    }
	    sm /= sw;
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		k = m[j + pp1 * m_dim1];
		z[j + z_dim1 * 3] -= sm * ty[k + js * ty_dim1];
/* L470: */
	    }
/* L480: */
	}
L490:
	sm = (float)0.;
	sv = sm;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    k = m[j + pp1 * m_dim1];
	    sm += w[k] * z[j + z_dim1 * 3];
	    z[k + (z_dim1 << 1)] = z[j + z_dim1];
/* L500: */
	}
	sm /= sw;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    z[j + z_dim1 * 3] -= sm;
/* Computing 2nd power */
	    d__1 = z[j + z_dim1 * 3];
	    sv += z[j + (z_dim1 << 2)] * (d__1 * d__1);
/* L510: */
	}
	sv /= sw;
	if (sv <= (float)0.) {
	    goto L520;
	}
	sv = (float)1. / sqrt(sv);
	goto L530;
L520:
	*ierr = 3;
	if (prams_1.itape > 0) {
	    io___26.ciunit = prams_1.itape;
	    s_wsfe(&io___26);
	    do_fio(&c__1, (char *)&is, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	return 0;
L530:
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    k = m[j + pp1 * m_dim1];
	    ty[k + is * ty_dim1] = z[j + z_dim1 * 3] * sv;
/* L540: */
	}
	sv = (float)0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
/* Computing 2nd power */
	    d__1 = ty[j + is * ty_dim1] - z[j + (z_dim1 << 1)];
	    sv += w[j] * (d__1 * d__1);
/* L550: */
	}
	rsq[is] = (float)1. - sv / sw;
	if (prams_1.itape > 0) {
	    io___27.ciunit = prams_1.itape;
	    s_wsfe(&io___27);
	    do_fio(&c__1, (char *)&iter, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&rsq[is], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	nt = nt % prams_1.nterm + 1;
	ct[nt - 1] = rsq[is];
	cmn = (float)100.;
	cmx = (float)-100.;
	i__2 = prams_1.nterm;
	for (i = 1; i <= i__2; ++i) {
/* Computing MIN */
	    d__1 = cmn, d__2 = ct[i - 1];
	    cmn = (real) min(d__1,d__2);
/* Computing MAX */
	    d__1 = cmx, d__2 = ct[i - 1];
	    cmx = (real) max(d__1,d__2);
/* L560: */
	}
	if (cmx - cmn <= *delrsq || iter >= prams_1.maxit) {
	    goto L570;
	}
	goto L330;
L570:
	if (prams_1.itape > 0) {
	    io___30.ciunit = prams_1.itape;
	    s_wsfe(&io___30);
	    do_fio(&c__1, (char *)&is, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&rsq[is], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
/* L580: */
    }
    return 0;
} /* mace_ */

/* Subroutine */ int model_(p, n, y, w, l, tx, ty, f, t, m, z)
integer *p, *n;
doublereal *y, *w;
integer *l;
doublereal *tx, *ty, *f, *t;
integer *m;
doublereal *z;
{
    /* System generated locals */
    integer m_dim1, m_offset, tx_dim1, tx_offset, z_dim1, z_offset, i__1, 
	    i__2;

    /* Local variables */
    extern /* Subroutine */ int sort_();
    static integer i, j, k;
    static doublereal s;
    static integer j1, j2;
    extern /* Subroutine */ int smothr_();
    static integer pp1;


/*          subroutine model(p,n,y,w,l,tx,ty,f,t,m,z) */
/* -------------------------------------------------------------------- */

/* computes response predictive  function f for the model yhat = f(t), */
/* where */
/*                                        p */
/*            f(t) = e(y : t),     t =   sum  tx<i> ( x<i> ) */
/*                                       i=1 */
/* using the x transformations tx constructed by subroutine ace. */
/* if y is a categorical variable (classification) then */
/*                                -1 */
/*                       f(t) = ty  (t). */
/* input: */

/*    p,n,y,w,l : same input as for subroutine ace. */
/*    tx,ty,m,z : output from subroutine ace. */

/* output: */

/*    f(n),t(n) : input for subroutine acemod. */

/* note: this subroutine must be called before subroutine acemod. */

/* ------------------------------------------------------------------- */

    /* Parameter adjustments */
    z_dim1 = *n;
    z_offset = z_dim1 + 1;
    z -= z_offset;
    m_dim1 = *n;
    m_offset = m_dim1 + 1;
    m -= m_offset;
    --t;
    --f;
    --ty;
    tx_dim1 = *n;
    tx_offset = tx_dim1 + 1;
    tx -= tx_offset;
    --w;
    --y;
    --l;

    /* Function Body */
    pp1 = *p + 1;
    if ((i__1 = l[pp1], abs(i__1)) != 5) {
	goto L20;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	t[j] = ty[j];
	m[j + pp1 * m_dim1] = j;
/* L10: */
    }
    goto L50;
L20:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	s = (float)0.;
	i__2 = *p;
	for (i = 1; i <= i__2; ++i) {
	    s += tx[j + i * tx_dim1];
/* L30: */
	}
	t[j] = s;
	m[j + pp1 * m_dim1] = j;
/* L40: */
    }
L50:
    sort_(&t[1], &m[pp1 * m_dim1 + 1], &c__1, n);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	k = m[j + pp1 * m_dim1];
	z[j + (z_dim1 << 1)] = w[k];
	if (y[k] >= prams_1.big) {
	    goto L60;
	}
	z[j + z_dim1] = y[k];
	goto L140;
L60:
	j1 = j;
	j2 = j1;
L70:
	if (y[m[j1 + pp1 * m_dim1]] < prams_1.big) {
	    goto L80;
	}
	--j1;
	if (j1 < 1) {
	    goto L80;
	}
	goto L70;
L80:
	if (y[m[j2 + pp1 * m_dim1]] < prams_1.big) {
	    goto L90;
	}
	++j2;
	if (j2 > *n) {
	    goto L90;
	}
	goto L80;
L90:
	if (j1 >= 1) {
	    goto L100;
	}
	k = j2;
	goto L130;
L100:
	if (j2 <= *n) {
	    goto L110;
	}
	k = j1;
	goto L130;
L110:
	if (t[j] - t[j1] >= t[j2] - t[j]) {
	    goto L120;
	}
	k = j1;
	goto L130;
L120:
	k = j2;
L130:
	z[j + z_dim1] = y[m[k + pp1 * m_dim1]];
	t[j] = t[k];
L140:
	;
    }
    if ((i__1 = l[pp1], abs(i__1)) != 5) {
	goto L160;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	f[j] = z[j + z_dim1];
/* L150: */
    }
    goto L170;
L160:
    smothr_(&c__1, n, &t[1], &z[z_offset], &z[(z_dim1 << 1) + 1], &f[1], &z[
	    z_dim1 * 6 + 1]);
L170:
    return 0;
} /* model_ */

/* Subroutine */ int acemod_(v, p, n, x, l, tx, f, t, m, yhat)
doublereal *v;
integer *p, *n;
doublereal *x;
integer *l;
doublereal *tx, *f, *t;
integer *m;
doublereal *yhat;
{
    /* System generated locals */
    integer m_dim1, m_offset, x_dim1, x_offset, tx_dim1, tx_offset, i__1, 
	    i__2;

    /* Local variables */
    static integer high, i, place, jh, jl;
    static doublereal th, vi, xt;
    static integer low;

/*          subroutine acemod(v,p,n,x,l,tx,f,t,m,yhat) */
/* -------------------------------------------------------------------- */

/* computes response y estimates from the model */

/*                yhat =  f ( t( v ) ) */

/* using the x transformations tx constructed by subroutine ace and */
/* the predictor function (f,t) constructed by subroutine model. */

/* input: */

/*       v(p) : vector of predictor values. */
/*    p,n,x,l : same input as for subroutine ace. */
/*       tx,m : output from subroutine ace. */
/*        f,t : output from subroutine model. */

/* output: */

/*    yhat : estimated response value for v. */

/* note: this subroutine must not be called before subroutine model. */

/* ------------------------------------------------------------------- */

    /* Parameter adjustments */
    --v;
    m_dim1 = *n;
    m_offset = m_dim1 + 1;
    m -= m_offset;
    --t;
    --f;
    tx_dim1 = *n;
    tx_offset = tx_dim1 + 1;
    tx -= tx_offset;
    x_dim1 = *p;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    --l;

    /* Function Body */
    th = (float)0.;
    i__1 = *p;
    for (i = 1; i <= i__1; ++i) {
	if (l[i] == 0) {
	    goto L90;
	}
	vi = v[i];
	if (vi < prams_1.big) {
	    goto L10;
	}
	if (x[i + m[*n + i * m_dim1] * x_dim1] >= prams_1.big) {
	    th += tx[m[*n + i * m_dim1] + i * tx_dim1];
	}
	goto L90;
L10:
	if (vi > x[i + m[i * m_dim1 + 1] * x_dim1]) {
	    goto L20;
	}
	place = 1;
	goto L80;
L20:
	if (vi < x[i + m[*n + i * m_dim1] * x_dim1]) {
	    goto L30;
	}
	place = *n;
	goto L80;
L30:
	low = 0;
	high = *n + 1;
L40:
	if (low + 1 >= high) {
	    goto L60;
	}
	place = (low + high) / 2;
	xt = x[i + m[place + i * m_dim1] * x_dim1];
	if (vi == xt) {
	    goto L80;
	}
	if (vi >= xt) {
	    goto L50;
	}
	high = place;
	goto L40;
L50:
	low = place;
	goto L40;
L60:
	if ((i__2 = l[i], abs(i__2)) == 5) {
	    goto L90;
	}
	jl = m[low + i * m_dim1];
	jh = m[high + i * m_dim1];
	if (x[i + jh * x_dim1] < prams_1.big) {
	    goto L70;
	}
	th += tx[jl + i * tx_dim1];
	goto L90;
L70:
	th = th + tx[jl + i * tx_dim1] + (tx[jh + i * tx_dim1] - tx[jl + i * 
		tx_dim1]) * (vi - x[i + jl * x_dim1]) / (x[i + jh * x_dim1] - 
		x[i + jl * x_dim1]);
	goto L90;
L80:
	th += tx[m[place + i * m_dim1] + i * tx_dim1];
L90:
	;
    }
    if (th > t[1]) {
	goto L100;
    }
    *yhat = f[1];
    return 0;
L100:
    if (th < t[*n]) {
	goto L110;
    }
    *yhat = f[*n];
    return 0;
L110:
    low = 0;
    high = *n + 1;
L120:
    if (low + 1 >= high) {
	goto L150;
    }
    place = (low + high) / 2;
    xt = t[place];
    if (th != xt) {
	goto L130;
    }
    *yhat = f[place];
    return 0;
L130:
    if (th >= xt) {
	goto L140;
    }
    high = place;
    goto L120;
L140:
    low = place;
    goto L120;
L150:
    if ((i__1 = l[*p + 1], abs(i__1)) != 5) {
	goto L170;
    }
    if (th - t[low] > t[high] - th) {
	goto L160;
    }
    *yhat = f[low];
    goto L180;
L160:
    *yhat = f[high];
    goto L180;
L170:
    *yhat = f[low] + (f[high] - f[low]) * (th - t[low]) / (t[high] - t[low]);
L180:
    return 0;
} /* acemod_ */


/*     block data */
/*     common /prams/ itape,maxit,nterm,span,alpha,big */

/* ------------------------------------------------------------------ */

/* these procedure parameters can be changed in the calling routine */
/* by defining the above labeled common and resetting the values with */
/* executable statements. */

/* itape : fortran file number for printer output. */
/*         (itape.le.0 => no printer output.) */
/* maxit : maximum number of iterations. */
/* nterm : number of consecutive iterations for which */
/*         rsq must change less than delcor for convergence. */
/* span, alpha : super smoother parameters (see below). */
/* big : a large representable floating point number. */

/* ------------------------------------------------------------------ */


/* Subroutine */ int scale_(p, n, w, sw, ty, tx, eps, maxit, r, sc)
integer *p, *n;
doublereal *w, *sw, *ty, *tx, *eps;
integer *maxit;
doublereal *r, *sc;
{
    /* System generated locals */
    integer tx_dim1, tx_offset, sc_dim1, sc_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static doublereal gama;
    static integer iter;
    static doublereal h;
    static integer i, j;
    static doublereal s, t, u, v, delta;
    static integer nit;

    /* Parameter adjustments */
    sc_dim1 = *p;
    sc_offset = sc_dim1 + 1;
    sc -= sc_offset;
    --r;
    tx_dim1 = *n;
    tx_offset = tx_dim1 + 1;
    tx -= tx_offset;
    --ty;
    --w;

    /* Function Body */
    i__1 = *p;
    for (i = 1; i <= i__1; ++i) {
	sc[i + sc_dim1] = (float)0.;
/* L10: */
    }
    nit = 0;
L20:
    ++nit;
    i__1 = *p;
    for (i = 1; i <= i__1; ++i) {
	sc[i + sc_dim1 * 5] = sc[i + sc_dim1];
/* L30: */
    }
    i__1 = *p;
    for (iter = 1; iter <= i__1; ++iter) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    s = (float)0.;
	    i__3 = *p;
	    for (i = 1; i <= i__3; ++i) {
		s += sc[i + sc_dim1] * tx[j + i * tx_dim1];
/* L40: */
	    }
	    r[j] = (ty[j] - s) * w[j];
/* L50: */
	}
	i__2 = *p;
	for (i = 1; i <= i__2; ++i) {
	    s = (float)0.;
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		s += r[j] * tx[j + i * tx_dim1];
/* L60: */
	    }
	    sc[i + (sc_dim1 << 1)] = s * (float)-2. / *sw;
/* L70: */
	}
	s = (float)0.;
	i__2 = *p;
	for (i = 1; i <= i__2; ++i) {
/* Computing 2nd power */
	    d__1 = sc[i + (sc_dim1 << 1)];
	    s += d__1 * d__1;
/* L80: */
	}
	if (s <= (float)0.) {
	    goto L170;
	}
	if (iter != 1) {
	    goto L100;
	}
	i__2 = *p;
	for (i = 1; i <= i__2; ++i) {
	    sc[i + sc_dim1 * 3] = -sc[i + (sc_dim1 << 1)];
/* L90: */
	}
	h = s;
	goto L120;
L100:
	gama = s / h;
	h = s;
	i__2 = *p;
	for (i = 1; i <= i__2; ++i) {
	    sc[i + sc_dim1 * 3] = -sc[i + (sc_dim1 << 1)] + gama * sc[i + (
		    sc_dim1 << 2)];
/* L110: */
	}
L120:
	s = (float)0.;
	t = s;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    u = (float)0.;
	    i__3 = *p;
	    for (i = 1; i <= i__3; ++i) {
		u += sc[i + sc_dim1 * 3] * tx[j + i * tx_dim1];
/* L130: */
	    }
	    s += u * r[j];
/* Computing 2nd power */
	    d__1 = u;
	    t += w[j] * (d__1 * d__1);
/* L140: */
	}
	delta = s / t;
	i__2 = *p;
	for (i = 1; i <= i__2; ++i) {
	    sc[i + sc_dim1] += delta * sc[i + sc_dim1 * 3];
	    sc[i + (sc_dim1 << 2)] = sc[i + sc_dim1 * 3];
/* L150: */
	}
/* L160: */
    }
L170:
    v = (float)0.;
    i__1 = *p;
    for (i = 1; i <= i__1; ++i) {
/* Computing MAX */
	d__2 = v, d__3 = (d__1 = sc[i + sc_dim1] - sc[i + sc_dim1 * 5], abs(
		d__1));
	v = (real) max(d__2,d__3);
/* L180: */
    }
    if (v < *eps || nit >= *maxit) {
	goto L190;
    }
    goto L20;
L190:
    i__1 = *p;
    for (i = 1; i <= i__1; ++i) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    tx[j + i * tx_dim1] = sc[i + sc_dim1] * tx[j + i * tx_dim1];
/* L200: */
	}
/* L210: */
    }
    return 0;
} /* scale_ */

