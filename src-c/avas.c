/* avas.f -- translated by f2c (version 19950110).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct parms_1_ {
    integer itape, maxit, nterm;
    real span, alpha;
};

#define parms_1 (*(struct parms_1_ *) &parms_)

struct spans_1_ {
    real spans[3];
};

#define spans_1 (*(struct spans_1_ *) &spans_)

struct consts_1_ {
    real big, sml, eps;
};

#define consts_1 (*(struct consts_1_ *) &consts_)

/* Initialized data */

struct {
    integer e_1[3];
    real e_2[2];
    } parms_ = { -6, 20, 3, (float)0., (float)5. };

struct {
    real e_1[3];
    } spans_ = { (float).05, (float).2, (float).5 };

struct {
    real e_1[3];
    } consts_ = { (float)1e20, (float)1e-4, (float).001 };


/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int avas_(p, n, x, y, w, l, delrsq, tx, ty, rsq, ierr, m, z, 
	yspan, iter, iters)
integer *p, *n;
doublereal *x, *y, *w;
integer *l;
doublereal *delrsq, *tx, *ty, *rsq;
integer *ierr, *m;
doublereal *z, *yspan;
integer *iter;
doublereal *iters;
{
    /* System generated locals */
    integer m_dim1, m_offset, x_dim1, x_offset, tx_dim1, tx_offset, z_dim1, 
	    z_offset, i__1, i__2;
    real r__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(), log(), exp();

    /* Local variables */
    static real rnew, tres;
    extern /* Subroutine */ int sort_();
    static integer i, j, k;
    extern /* Subroutine */ int ctsub_(), rlsmo_();
    static doublereal ct[10];
    static integer np;
    static doublereal sm;
    extern /* Subroutine */ int bakfit_();
    static integer nt;
    extern /* Subroutine */ int calcmu_();
    static real rr;
    static doublereal sv, sw;
    static real sumlog;
    static integer pp1, pp2;
    static real dof, cmn, cmx, rss;
    static doublereal svx;

    /* Parameter adjustments */
    z_dim1 = *n;
    z_offset = z_dim1 + 1;
    z -= z_offset;
    m_dim1 = *n;
    m_offset = m_dim1 + 1;
    m -= m_offset;
    --ty;
    tx_dim1 = *n;
    tx_offset = tx_dim1 + 1;
    tx -= tx_offset;
    --w;
    --y;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    --l;
    iters -= 101;

    /* Function Body */
    *ierr = 0;
    pp1 = *p + 1;
    pp2 = *p + 2;
    sm = (float)0.;
    sv = sm;
    sw = sv;
    np = 0;
    i__1 = *p;
    for (i = 1; i <= i__1; ++i) {
	if (! (l[i] > 0)) {
	    goto L23002;
	}
	++np;
L23002:
/* L23000: */
	;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	sm += w[j] * y[j];
/* Computing 2nd power */
	d__1 = y[j];
	sv += w[j] * (d__1 * d__1);
	sw += w[j];
	m[j + pp1 * m_dim1] = j;
	z[j + (z_dim1 << 1)] = y[j];
/* L23004: */
    }
    sm /= sw;
/* Computing 2nd power */
    d__1 = sm;
    sv = sv / sw - d__1 * d__1;
    sv = (float)1. / sqrt(sv);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z[j + z_dim1] = (y[j] - sm) * sv;
/* L23006: */
    }
    sort_(&z[(z_dim1 << 1) + 1], &m[pp1 * m_dim1 + 1], &c__1, n);
    i__1 = *p;
    for (i = 1; i <= i__1; ++i) {
	if (! (l[i] > 0)) {
	    goto L23010;
	}
	sm = 0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    sm += w[j] * x[j + i * x_dim1];
/* L23012: */
	}
	sm /= sw;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    m[j + i * m_dim1] = j;
	    z[j + (z_dim1 << 1)] = x[j + i * x_dim1];
/* L23014: */
	}
	sort_(&z[(z_dim1 << 1) + 1], &m[i * m_dim1 + 1], &c__1, n);
L23010:
/* L23008: */
	;
    }
    *rsq = (float)0.;
    *iter = 0;
    parms_1.nterm = min(parms_1.nterm,10);
    nt = 0;
    i__1 = parms_1.nterm;
    for (i = 1; i <= i__1; ++i) {
	ct[i - 1] = (float)100.;
/* L23016: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	ty[j] = z[j + z_dim1];
/* L23018: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z[j + z_dim1 * 9] = ty[j];
/* L23020: */
    }
    bakfit_(iter, delrsq, rsq, &sw, &l[1], &z[z_offset], &m[m_offset], &x[
	    x_offset], &z[z_dim1 * 9 + 1], &tx[tx_offset], &w[1], n, p, &np);
    sumlog = (float)0.;
L23022:
    ++(*iter);
    if (! (l[pp1] == 4)) {
	goto L23025;
    }
    goto L992;
L23025:
    calcmu_(n, p, &l[1], &z[z_offset], &tx[tx_offset]);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	tres = ty[j] - z[j + z_dim1 * 10];
	if (! (dabs(tres) < (float)1e-10)) {
	    goto L23029;
	}
	tres = (float)1e-10;
L23029:
/* Computing 2nd power */
	r__1 = tres;
	z[j + (z_dim1 << 1)] = log(sqrt(r__1 * r__1));
	m[j + pp2 * m_dim1] = j;
/* L23027: */
    }
    sort_(&z[z_dim1 * 10 + 1], &m[pp2 * m_dim1 + 1], &c__1, n);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	k = m[j + pp2 * m_dim1];
	z[j + (z_dim1 << 2)] = z[k + (z_dim1 << 1)];
	z[j + z_dim1 * 5] = w[k];
/* L23031: */
    }
    rlsmo_(&z[z_dim1 * 10 + 1], &z[(z_dim1 << 2) + 1], &z[z_dim1 * 5 + 1], 
	    yspan, &dof, n, &z[z_dim1 * 6 + 1], &rss, &z[z_dim1 * 7 + 1]);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	k = m[j + pp2 * m_dim1];
	z[j + z_dim1 * 7] = exp(-z[j + z_dim1 * 6]);
	sumlog += *n * (w[j] / sw) * 2 * z[j + z_dim1 * 6];
	z[j + (z_dim1 << 3)] = ty[k];
/* L23033: */
    }
    ctsub_(n, &z[z_dim1 * 10 + 1], &z[z_dim1 * 7 + 1], &z[(z_dim1 << 3) + 1], 
	    &z[z_dim1 * 9 + 1]);
    sm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	sm += w[j] * z[j + z_dim1 * 9];
/* L23035: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	k = m[j + pp2 * m_dim1];
	ty[k] = z[j + z_dim1 * 9] - sm / sw;
/* L23037: */
    }
    sv = 0.;
    svx = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	sv += w[j] / sw * ty[j] * ty[j];
	svx += w[j] / sw * z[j + z_dim1 * 10] * z[j + z_dim1 * 10];
/* L23039: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	ty[j] /= sqrt(sv);
	i__2 = *p;
	for (i = 1; i <= i__2; ++i) {
	    if (! (l[i] > 0)) {
		goto L23045;
	    }
	    tx[j + i * tx_dim1] /= sqrt(svx);
L23045:
/* L23043: */
	    ;
	}
/* L23041: */
    }
L992:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z[j + z_dim1 * 9] = ty[j];
/* L23047: */
    }
    bakfit_(iter, delrsq, rsq, &sw, &l[1], &z[z_offset], &m[m_offset], &x[
	    x_offset], &z[z_dim1 * 9 + 1], &tx[tx_offset], &w[1], n, p, &np);
    sumlog += *n * log(sv);
    rr = (float)0.;
    calcmu_(n, p, &l[1], &z[z_offset], &tx[tx_offset]);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
	d__1 = ty[j] - z[j + z_dim1 * 10];
	rr += w[j] / sw * (d__1 * d__1);
/* L23049: */
    }
    *rsq = 1 - rr;
    rnew = sumlog + rr;
    iters[*iter + 100] = (doublereal) (*iter);
    iters[*iter + 200] = *rsq;
    nt = nt % parms_1.nterm + 1;
    ct[nt - 1] = *rsq;
    cmn = (float)100.;
    cmx = (float)-100.;
    i__1 = parms_1.nterm;
    for (i = 1; i <= i__1; ++i) {
/* Computing MIN */
	d__1 = cmn, d__2 = ct[i - 1];
	cmn = (real) min(d__1,d__2);
/* Computing MAX */
	d__1 = cmx, d__2 = ct[i - 1];
	cmx = (real) max(d__1,d__2);
/* L23051: */
    }
    if (! (cmx - cmn <= *delrsq || *iter >= parms_1.maxit || l[pp1] == 4)) {
	goto L23053;
    }
    return 0;
L23053:
/* L23023: */
    goto L23022;
    return 0;
} /* avas_ */

/* Subroutine */ int calcmu_(n, p, l, z, tx)
integer *n, *p, *l;
doublereal *z, *tx;
{
    /* System generated locals */
    integer z_dim1, z_offset, tx_dim1, tx_offset, i__1, i__2;

    /* Local variables */
    static integer i, j;

    /* Parameter adjustments */
    z_dim1 = *n;
    z_offset = z_dim1 + 1;
    z -= z_offset;
    tx_dim1 = *n;
    tx_offset = tx_dim1 + 1;
    tx -= tx_offset;
    --l;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z[j + z_dim1 * 10] = 0.;
	i__2 = *p;
	for (i = 1; i <= i__2; ++i) {
	    if (! (l[i] > 0)) {
		goto L23059;
	    }
	    z[j + z_dim1 * 10] += tx[j + i * tx_dim1];
L23059:
/* L23057: */
	    ;
	}
/* L23055: */
    }
    return 0;
} /* calcmu_ */

/* Subroutine */ int bakfit_(iter, delrsq, rsq, sw, l, z, m, x, ty, tx, w, n, 
	p, np)
integer *iter;
doublereal *delrsq, *rsq, *sw;
integer *l;
doublereal *z;
integer *m;
doublereal *x, *ty, *tx, *w;
integer *n, *p, *np;
{
    /* System generated locals */
    integer m_dim1, m_offset, z_dim1, z_offset, tx_dim1, tx_offset, x_dim1, 
	    x_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static real rsqi;
    static integer i, j, k;
    static doublereal sm;
    extern /* Subroutine */ int calcmu_();
    static doublereal sv;
    extern /* Subroutine */ int smothr_();
    static integer nit;

    /* Parameter adjustments */
    --l;
    --w;
    --ty;
    m_dim1 = *n;
    m_offset = m_dim1 + 1;
    m -= m_offset;
    z_dim1 = *n;
    z_offset = z_dim1 + 1;
    z -= z_offset;
    tx_dim1 = *n;
    tx_offset = tx_dim1 + 1;
    tx -= tx_offset;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;

    /* Function Body */
    calcmu_(n, p, &l[1], &z[z_offset], &tx[tx_offset]);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	ty[j] -= z[j + z_dim1 * 10];
/* L23061: */
    }
    nit = 0;
L23063:
    rsqi = *rsq;
    ++nit;
    i__1 = *p;
    for (i = 1; i <= i__1; ++i) {
	if (! (l[i] > 0)) {
	    goto L23068;
	}
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    k = m[j + i * m_dim1];
	    z[j + z_dim1] = ty[k] + tx[k + i * tx_dim1];
	    z[j + (z_dim1 << 1)] = x[k + i * x_dim1];
	    z[j + z_dim1 * 7] = w[k];
/* L23070: */
	}
	smothr_(&l[i], n, &z[(z_dim1 << 1) + 1], &z[z_offset], &z[z_dim1 * 7 
		+ 1], &z[z_dim1 * 6 + 1], &z[z_dim1 * 11 + 1]);
	sm = (float)0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    sm += z[j + z_dim1 * 7] * z[j + z_dim1 * 6];
/* L23072: */
	}
	sm /= *sw;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    z[j + z_dim1 * 6] -= sm;
/* L23074: */
	}
	sv = (float)0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
/* Computing 2nd power */
	    d__1 = z[j + z_dim1] - z[j + z_dim1 * 6];
	    sv += z[j + z_dim1 * 7] * (d__1 * d__1);
/* L23076: */
	}
	sv = (float)1. - sv / *sw;
	*rsq = sv;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    k = m[j + i * m_dim1];
	    tx[k + i * tx_dim1] = z[j + z_dim1 * 6];
	    ty[k] = z[j + z_dim1] - z[j + z_dim1 * 6];
/* L23078: */
	}
L23068:
/* L23066: */
	;
    }
/* L23064: */
    if (! (*np == 1 || (d__1 = *rsq - rsqi, abs(d__1)) <= *delrsq || nit >= 
	    parms_1.maxit)) {
	goto L23063;
    }
    if (! (*rsq == (float)0. && *iter == 0)) {
	goto L23080;
    }
    i__1 = *p;
    for (i = 1; i <= i__1; ++i) {
	if (! (l[i] > 0)) {
	    goto L23084;
	}
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    tx[j + i * tx_dim1] = x[j + i * x_dim1];
/* L23086: */
	}
L23084:
/* L23082: */
	;
    }
L23080:
    return 0;
} /* bakfit_ */

/* Subroutine */ int ctsub_(n, u, v, y, ty)
integer *n;
doublereal *u, *v, *y, *ty;
{
    static integer i, j;

    /* Parameter adjustments */
    --ty;
    --y;
    --v;
    --u;

    /* Function Body */
    i = 1;
L23088:
    if (! (i <= *n)) {
	goto L23090;
    }
    if (! (y[i] <= u[1])) {
	goto L23091;
    }
    ty[i] = (y[i] - u[1]) * v[1];
    goto L23092;
L23091:
    j = 1;
    ty[i] = 0.;
L23093:
    if (! (j <= *n && y[i] > u[j])) {
	goto L23094;
    }
    if (! (j > 1)) {
	goto L23095;
    }
    ty[i] += (u[j] - u[j - 1]) * (v[j] + v[j - 1]) / 2;
L23095:
    ++j;
    goto L23093;
L23094:
    if (! (y[i] <= u[*n])) {
	goto L23097;
    }
    ty[i] += (y[i] - u[j - 1]) * (float).5 * (v[j - 1] * 2 + (y[i] - u[j - 1])
	     * (v[j] - v[j - 1]) / (u[j] - u[j - 1]));
    goto L23098;
L23097:
    ty[i] += (y[i] - u[*n]) * v[*n];
L23098:
L23092:
    ++i;
    goto L23088;
L23090:
    return 0;
} /* ctsub_ */

/* ------------------------------------------------------------------ */

/* these procedure parameters can be changed in the calling routine */
/* by defining the above labeled common and resetting the values with */
/* executable statements. */

/* itape : fortran file number for printer output. */
/*         (itape.le.0 => no printer output.) */
/* maxit : maximum number of iterations. */
/* nterm : number of consecutive iterations for which */
/*         rsq must change less than delcor for convergence. */
/* span, alpha : super smoother parameters. */
/*   (see - friedman and stuetzle, reference above.) */

/* ------------------------------------------------------------------ */
/* --------------------------------------------------------------- */

/* this sets the compile time (default) values for various */
/* internal parameters : */

/* spans : span values for the three running linear smoothers. */
/* spans(1) : tweeter span. */
/* spans(2) : midrange span. */
/* spans(3) : woofer span. */
/* (these span values should be changed only with care.) */
/* big : a large representable floating point number. */
/* sml : a small number. should be set so that (sml)**(10.0) does */
/*       not cause floating point underflow. */
/* eps : used to numerically stabilize slope calculations for */
/*       running linear fits. */

/* these parameter values can be changed by declaring the */
/* relevant labeled common in the main program and resetting */
/* them with executable statements. */

/* ----------------------------------------------------------------- */

/* Subroutine */ int smothr_(l, n, x, y, w, smo, scr)
integer *l, *n;
doublereal *x, *y, *w, *smo, *scr;
{
    /* System generated locals */
    integer scr_dim1, scr_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static doublereal a, b, d;
    static integer i, j, j0;
    static doublereal sm, sw;
    extern /* Subroutine */ int montne_(), supsmu_();

    /* Parameter adjustments */
    scr_dim1 = *n;
    scr_offset = scr_dim1 + 1;
    scr -= scr_offset;
    --smo;
    --w;
    --y;
    --x;

    /* Function Body */
    if (*l < 5) {
	goto L50;
    }
    j = 1;
L10:
    j0 = j;
    sm = w[j] * y[j];
    sw = w[j];
    if (j >= *n) {
	goto L30;
    }
L20:
    if (x[j + 1] > x[j]) {
	goto L30;
    }
    ++j;
    sm += w[j] * y[j];
    sw += w[j];
    if (j >= *n) {
	goto L30;
    }
    goto L20;
L30:
    sm /= sw;
    i__1 = j;
    for (i = j0; i <= i__1; ++i) {
	smo[i] = sm;
/* L40: */
    }
    ++j;
    if (j > *n) {
	goto L250;
    }
    goto L10;
L50:
    if (*l != 4) {
	goto L80;
    }
    sm = (float)0.;
    sw = sm;
    b = sw;
    d = b;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	sm += w[j] * x[j] * y[j];
/* Computing 2nd power */
	d__1 = x[j];
	sw += w[j] * (d__1 * d__1);
	b += w[j] * x[j];
	d += w[j];
/* L60: */
    }
/* Computing 2nd power */
    d__1 = b;
    a = sm / (sw - d__1 * d__1 / d);
    b /= d;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	smo[j] = a * (x[j] - b);
/* L70: */
    }
    goto L250;
L80:
    supsmu_(n, &x[1], &y[1], &w[1], l, &parms_1.span, &parms_1.alpha, &smo[1],
	     &scr[scr_offset]);
    if (*l != 3) {
	goto L250;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	scr[j + scr_dim1] = smo[j];
	scr[*n - j + 1 + (scr_dim1 << 1)] = scr[j + scr_dim1];
/* L90: */
    }
    montne_(&scr[scr_offset], n);
    montne_(&scr[(scr_dim1 << 1) + 1], n);
    sm = (float)0.;
    sw = sm;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
	d__1 = smo[j] - scr[j + scr_dim1];
	sm += d__1 * d__1;
/* Computing 2nd power */
	d__1 = smo[j] - scr[*n - j + 1 + (scr_dim1 << 1)];
	sw += d__1 * d__1;
/* L100: */
    }
    if (sm >= sw) {
	goto L120;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	smo[j] = scr[j + scr_dim1];
/* L110: */
    }
    goto L140;
L120:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	smo[j] = scr[*n - j + 1 + (scr_dim1 << 1)];
/* L130: */
    }
L140:
    j = 1;
L150:
    j0 = j;
    if (j >= *n) {
	goto L170;
    }
L160:
    if (smo[j + 1] != smo[j]) {
	goto L170;
    }
    ++j;
    if (j >= *n) {
	goto L170;
    }
    goto L160;
L170:
    if (j <= j0) {
	goto L190;
    }
    a = (float)0.;
    if (j0 > 1) {
	a = (smo[j0] - smo[j0 - 1]) * (float).5;
    }
    b = (float)0.;
    if (j < *n) {
	b = (smo[j + 1] - smo[j]) * (float).5;
    }
    d = (a + b) / (j - j0);
    if (a == (float)0. || b == (float)0.) {
	d *= (float)2.;
    }
    if (a == (float)0.) {
	a = b;
    }
    i__1 = j;
    for (i = j0; i <= i__1; ++i) {
	smo[i] = smo[i] - a + d * (i - j0);
/* L180: */
    }
L190:
    ++j;
    if (j > *n) {
	goto L200;
    }
    goto L150;
L200:
    j = 1;
L210:
    j0 = j;
    sm = smo[j];
    if (j >= *n) {
	goto L230;
    }
L220:
    if (x[j + 1] > x[j]) {
	goto L230;
    }
    ++j;
    sm += smo[j];
    if (j >= *n) {
	goto L230;
    }
    goto L220;
L230:
    sm /= j - j0 + 1;
    i__1 = j;
    for (i = j0; i <= i__1; ++i) {
	smo[i] = sm;
/* L240: */
    }
    ++j;
    if (j > *n) {
	goto L250;
    }
    goto L210;
L250:
    return 0;
} /* smothr_ */

/* Subroutine */ int montne_(x, n)
doublereal *x;
integer *n;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, bb, eb, bl, el, br, er;
    static real pmn;

    /* Parameter adjustments */
    --x;

    /* Function Body */
    bb = 0;
    eb = bb;
L10:
    if (eb >= *n) {
	goto L110;
    }
    bb = eb + 1;
    eb = bb;
L20:
    if (eb >= *n) {
	goto L30;
    }
    if (x[bb] != x[eb + 1]) {
	goto L30;
    }
    ++eb;
    goto L20;
L30:
    if (eb >= *n) {
	goto L70;
    }
    if (x[eb] <= x[eb + 1]) {
	goto L70;
    }
    br = eb + 1;
    er = br;
L40:
    if (er >= *n) {
	goto L50;
    }
    if (x[er + 1] != x[br]) {
	goto L50;
    }
    ++er;
    goto L40;
L50:
    pmn = (x[bb] * (eb - bb + 1) + x[br] * (er - br + 1)) / (er - bb + 1);
    eb = er;
    i__1 = eb;
    for (i = bb; i <= i__1; ++i) {
	x[i] = pmn;
/* L60: */
    }
L70:
    if (bb <= 1) {
	goto L10;
    }
    if (x[bb - 1] <= x[bb]) {
	goto L10;
    }
    bl = bb - 1;
    el = bl;
L80:
    if (bl <= 1) {
	goto L90;
    }
    if (x[bl - 1] != x[el]) {
	goto L90;
    }
    --bl;
    goto L80;
L90:
    pmn = (x[bb] * (eb - bb + 1) + x[bl] * (el - bl + 1)) / (eb - bl + 1);
    bb = bl;
    i__1 = eb;
    for (i = bb; i <= i__1; ++i) {
	x[i] = pmn;
/* L100: */
    }
    goto L30;
L110:
    return 0;
} /* montne_ */

/* Subroutine */ int sort_(v, a, ii, jj)
doublereal *v;
integer *a, *ii, *jj;
{
    static integer i, j, k, l, m, t, ij, il[20], iu[20], tt;
    static real vt, vtt;


/*     puts into a the permutation vector which sorts v into */
/*     increasing order.  only elements from ii to jj are considered. */
/*     arrays iu(k) and il(k) permit sorting up to 2**(k+1)-1 elements */

/*     this is a modification of cacm algorithm #347 by r. c. singleton, 
*/
/*     which is a modified hoare quicksort. */

    /* Parameter adjustments */
    --v;
    --a;

    /* Function Body */
    m = 1;
    i = *ii;
    j = *jj;
L10:
    if (i >= j) {
	goto L80;
    }
L20:
    k = i;
    ij = (j + i) / 2;
    t = a[ij];
    vt = v[ij];
    if (v[i] <= vt) {
	goto L30;
    }
    a[ij] = a[i];
    a[i] = t;
    t = a[ij];
    v[ij] = v[i];
    v[i] = vt;
    vt = v[ij];
L30:
    l = j;
    if (v[j] >= vt) {
	goto L50;
    }
    a[ij] = a[j];
    a[j] = t;
    t = a[ij];
    v[ij] = v[j];
    v[j] = vt;
    vt = v[ij];
    if (v[i] <= vt) {
	goto L50;
    }
    a[ij] = a[i];
    a[i] = t;
    t = a[ij];
    v[ij] = v[i];
    v[i] = vt;
    vt = v[ij];
    goto L50;
L40:
    a[l] = a[k];
    a[k] = tt;
    v[l] = v[k];
    v[k] = vtt;
L50:
    --l;
    if (v[l] > vt) {
	goto L50;
    }
    tt = a[l];
    vtt = v[l];
L60:
    ++k;
    if (v[k] < vt) {
	goto L60;
    }
    if (k <= l) {
	goto L40;
    }
    if (l - i <= j - k) {
	goto L70;
    }
    il[m - 1] = i;
    iu[m - 1] = l;
    i = k;
    ++m;
    goto L90;
L70:
    il[m - 1] = k;
    iu[m - 1] = j;
    j = l;
    ++m;
    goto L90;
L80:
    --m;
    if (m == 0) {
	return 0;
    }
    i = il[m - 1];
    j = iu[m - 1];
L90:
    if (j - i > 10) {
	goto L20;
    }
    if (i == *ii) {
	goto L10;
    }
    --i;
L100:
    ++i;
    if (i == j) {
	goto L80;
    }
    t = a[i + 1];
    vt = v[i + 1];
    if (v[i] <= vt) {
	goto L100;
    }
    k = i;
L110:
    a[k + 1] = a[k];
    v[k + 1] = v[k];
    --k;
    if (vt < v[k]) {
	goto L110;
    }
    a[k + 1] = t;
    v[k + 1] = vt;
    goto L100;
} /* sort_ */

/* Subroutine */ int supsmu_(n, x, y, w, iper, span, alpha, smo, sc)
integer *n;
doublereal *x, *y, *w;
integer *iper;
real *span, *alpha;
doublereal *smo, *sc;
{
    /* System generated locals */
    integer sc_dim1, sc_offset, i__1;
    real r__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double pow_dd();

    /* Local variables */
    static integer jper;
    static real a, f;
    static doublereal h;
    static integer i, j;
    static real scale, sw, sy, resmin;
    extern /* Subroutine */ int smooth_();
    static real vsmlsq;

/* ------------------------------------------------------------------ */

/* super smoother (friedman and stuetzle, 1984). */

/* version 3/10/84 */

/* coded by: j. h. friedman */
/*           department of statistics and */
/*           stanford linear accelerator center */
/*           stanford university */
/*           stanford ca. 94305 */

/* input: */
/*    n : number of observations (x,y - pairs). */
/*    x(n) : ordered abscissa values. */
/*    y(n) : corresponding ordinate (response) values. */
/*    w(n) : weight for each (x,y) observation. */
/*    iper : periodic variable flag. */
/*       iper=1 => x is ordered interval variable. */
/*       iper=2 => x is a periodic variable with values */
/*                 in the range (0.0,1.0) and peroid 1.0. */
/*    span : smoother span (fraction of observations in window). */
/*           span=0.0 => automatic (variable) span selection. */
/*    alpha : controles high frequency (small span) penality */
/*            used with automatic span selection (base tone control). */
/*            (alpha.le.0.0 or alpha.gt.10.0 => no effect.) */
/* output: */
/*   smo(n) : smoothed ordinate (response) values. */
/* scratch: */
/*   sc(n,7) : internal working storage. */

/* note: */
/*    for small samples (n < 40) or if there are substantial serial */
/*    correlations between obserations close in x - value, then */
/*    a prespecified fixed span smoother (span > 0) should be */
/*    used. reasonable span values are 0.3 to 0.5. */

/* ------------------------------------------------------------------ */
    /* Parameter adjustments */
    sc_dim1 = *n;
    sc_offset = sc_dim1 + 1;
    sc -= sc_offset;
    --smo;
    --w;
    --y;
    --x;

    /* Function Body */
    if (x[*n] > x[1]) {
	goto L30;
    }
    sy = (float)0.;
    sw = sy;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	sy += w[j] * y[j];
	sw += w[j];
/* L10: */
    }
    a = sy / sw;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	smo[j] = a;
/* L20: */
    }
    return 0;
L30:
    i = *n / 4;
    j = i * 3;
    scale = x[j] - x[i];
L40:
    if (scale > (float)0.) {
	goto L50;
    }
    if (j < *n) {
	++j;
    }
    if (i > 1) {
	--i;
    }
    scale = x[j] - x[i];
    goto L40;
L50:
/* Computing 2nd power */
    r__1 = consts_1.eps * scale;
    vsmlsq = r__1 * r__1;
    jper = *iper;
    if (*iper == 2 && (x[1] < (float)0. || x[*n] > (float)1.)) {
	jper = 1;
    }
    if (jper < 1 || jper > 2) {
	jper = 1;
    }
    if (*span <= (float)0.) {
	goto L60;
    }
    smooth_(n, &x[1], &y[1], &w[1], span, &jper, &vsmlsq, &smo[1], &sc[
	    sc_offset]);
    return 0;
L60:
    for (i = 1; i <= 3; ++i) {
	smooth_(n, &x[1], &y[1], &w[1], &spans_1.spans[i - 1], &jper, &vsmlsq,
		 &sc[((i << 1) - 1) * sc_dim1 + 1], &sc[sc_dim1 * 7 + 1]);
	i__1 = -jper;
	smooth_(n, &x[1], &sc[sc_dim1 * 7 + 1], &w[1], &spans_1.spans[1], &
		i__1, &vsmlsq, &sc[(i << 1) * sc_dim1 + 1], &h);
/* L70: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	resmin = consts_1.big;
	for (i = 1; i <= 3; ++i) {
	    if (sc[j + (i << 1) * sc_dim1] >= resmin) {
		goto L80;
	    }
	    resmin = sc[j + (i << 1) * sc_dim1];
	    sc[j + sc_dim1 * 7] = spans_1.spans[i - 1];
L80:
	    ;
	}
	if (*alpha > (float)0. && *alpha <= (float)10. && resmin < sc[j + 
		sc_dim1 * 6]) {
/* Computing MAX */
	    d__2 = consts_1.sml, d__3 = resmin / sc[j + sc_dim1 * 6];
	    d__1 = (doublereal) ((real) max(d__2,d__3));
	    d__4 = (doublereal) ((float)10. - *alpha);
	    sc[j + sc_dim1 * 7] += (spans_1.spans[2] - sc[j + sc_dim1 * 7]) * 
		    pow_dd(&d__1, &d__4);
	}
/* L90: */
    }
    i__1 = -jper;
    smooth_(n, &x[1], &sc[sc_dim1 * 7 + 1], &w[1], &spans_1.spans[1], &i__1, &
	    vsmlsq, &sc[(sc_dim1 << 1) + 1], &h);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (sc[j + (sc_dim1 << 1)] <= spans_1.spans[0]) {
	    sc[j + (sc_dim1 << 1)] = spans_1.spans[0];
	}
	if (sc[j + (sc_dim1 << 1)] >= spans_1.spans[2]) {
	    sc[j + (sc_dim1 << 1)] = spans_1.spans[2];
	}
	f = sc[j + (sc_dim1 << 1)] - spans_1.spans[1];
	if (f >= (float)0.) {
	    goto L100;
	}
	f = -(doublereal)f / (spans_1.spans[1] - spans_1.spans[0]);
	sc[j + (sc_dim1 << 2)] = ((float)1. - f) * sc[j + sc_dim1 * 3] + f * 
		sc[j + sc_dim1];
	goto L110;
L100:
	f /= spans_1.spans[2] - spans_1.spans[1];
	sc[j + (sc_dim1 << 2)] = ((float)1. - f) * sc[j + sc_dim1 * 3] + f * 
		sc[j + sc_dim1 * 5];
L110:
	;
    }
    i__1 = -jper;
    smooth_(n, &x[1], &sc[(sc_dim1 << 2) + 1], &w[1], spans_1.spans, &i__1, &
	    vsmlsq, &smo[1], &h);
    return 0;
} /* supsmu_ */

/* Subroutine */ int smooth_(n, x, y, w, span, iper, vsmlsq, smo, acvr)
integer *n;
doublereal *x, *y, *w;
real *span;
integer *iper;
real *vsmlsq;
doublereal *smo, *acvr;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static real cvar;
    static integer jper;
    static real a, h;
    static integer i, j, j0, in, it;
    static real xm, ym, wt, sy, fbo, fbw;
    static integer ibw;
    static real var, tmp, xti;
    static integer out;
    static real xto;

    /* Parameter adjustments */
    --acvr;
    --smo;
    --w;
    --y;
    --x;

    /* Function Body */
    xm = (float)0.;
    ym = xm;
    var = ym;
    cvar = var;
    fbw = cvar;
    jper = abs(*iper);
    ibw = *span * (float).5 * *n + (float).5;
    if (ibw < 2) {
	ibw = 2;
    }
    it = (ibw << 1) + 1;
    i__1 = it;
    for (i = 1; i <= i__1; ++i) {
	j = i;
	if (jper == 2) {
	    j = i - ibw - 1;
	}
	xti = x[j];
	if (j >= 1) {
	    goto L10;
	}
	j = *n + j;
	xti = x[j] - (float)1.;
L10:
	wt = w[j];
	fbo = fbw;
	fbw += wt;
	xm = (fbo * xm + wt * xti) / fbw;
	ym = (fbo * ym + wt * y[j]) / fbw;
	tmp = (float)0.;
	if (fbo > (float)0.) {
	    tmp = fbw * wt * (xti - xm) / fbo;
	}
	var += tmp * (xti - xm);
	cvar += tmp * (y[j] - ym);
/* L20: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	out = j - ibw - 1;
	in = j + ibw;
	if (jper != 2 && (out < 1 || in > *n)) {
	    goto L60;
	}
	if (out >= 1) {
	    goto L30;
	}
	out = *n + out;
	xto = x[out] - (float)1.;
	xti = x[in];
	goto L50;
L30:
	if (in <= *n) {
	    goto L40;
	}
	in -= *n;
	xti = x[in] + (float)1.;
	xto = x[out];
	goto L50;
L40:
	xto = x[out];
	xti = x[in];
L50:
	wt = w[out];
	fbo = fbw;
	fbw -= wt;
	tmp = (float)0.;
	if (fbw > (float)0.) {
	    tmp = fbo * wt * (xto - xm) / fbw;
	}
	var -= tmp * (xto - xm);
	cvar -= tmp * (y[out] - ym);
	xm = (fbo * xm - wt * xto) / fbw;
	ym = (fbo * ym - wt * y[out]) / fbw;
	wt = w[in];
	fbo = fbw;
	fbw += wt;
	xm = (fbo * xm + wt * xti) / fbw;
	ym = (fbo * ym + wt * y[in]) / fbw;
	tmp = (float)0.;
	if (fbo > (float)0.) {
	    tmp = fbw * wt * (xti - xm) / fbo;
	}
	var += tmp * (xti - xm);
	cvar += tmp * (y[in] - ym);
L60:
	a = (float)0.;
	if (var > *vsmlsq) {
	    a = cvar / var;
	}
	smo[j] = a * (x[j] - xm) + ym;
	if (*iper <= 0) {
	    goto L70;
	}
	h = (float)1. / fbw;
	if (var > *vsmlsq) {
/* Computing 2nd power */
	    d__1 = x[j] - xm;
	    h += d__1 * d__1 / var;
	}
	acvr[j] = (d__1 = y[j] - smo[j], abs(d__1)) / ((float)1. - w[j] * h);
L70:
	;
    }
    j = 1;
L80:
    j0 = j;
    sy = smo[j] * w[j];
    fbw = w[j];
    if (j >= *n) {
	goto L100;
    }
L90:
    if (x[j + 1] > x[j]) {
	goto L100;
    }
    ++j;
    sy += w[j] * smo[j];
    fbw += w[j];
    if (j >= *n) {
	goto L100;
    }
    goto L90;
L100:
    if (j <= j0) {
	goto L120;
    }
    sy /= fbw;
    i__1 = j;
    for (i = j0; i <= i__1; ++i) {
	smo[i] = sy;
/* L110: */
    }
L120:
    ++j;
    if (j > *n) {
	goto L130;
    }
    goto L80;
L130:
    return 0;
} /* smooth_ */

