/* rlsmo.f -- translated by f2c (version 19950110).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;

/*     MORTRAN 2.79 (RESERVED KEYWORD MACROS OF 09/28/81) */
/* Subroutine */ int rlsmo_(x, y, w, span, dof, n, smo, rss, scrat)
doublereal *x, *y, *w, *span;
real *dof;
integer *n;
doublereal *smo, *rss, *scrat;
{
    /* Initialized data */

    static doublereal cvspan[6] = { .3,.4,.5,.6,.7,1. };

    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int smth_();
    static integer i, k;
    static real penal;
    static integer idmin;
    static doublereal cvmin;
    static integer cross;
    static doublereal cvrss[6];
    static real s0;

    /* Parameter adjustments */
    --scrat;
    --smo;
    --w;
    --y;
    --x;

    /* Function Body */
    cross = 0;
    if (*span == 0.) {
	cross = 1;
    }
    penal = (float).01;
    cvmin = (float)1e15;
    idmin = 1;
    if (cross != 1) {
	goto L10021;
    }
    k = 1;
    goto L10033;
L10031:
    ++k;
L10033:
    if (k > 6) {
	goto L10032;
    }
    smth_(&x[1], &y[1], &w[1], &cvspan[k - 1], dof, n, &c__1, &smo[1], &s0, &
	    cvrss[k - 1], &scrat[1]);
    if (cvrss[k - 1] > cvmin) {
	goto L10051;
    }
    cvmin = cvrss[k - 1];
    idmin = k;
L10051:
    goto L10031;
L10032:
    *span = cvspan[idmin - 1];
    if (penal <= (float)0.) {
	goto L10071;
    }
    cvmin = (penal + (float)1.) * cvmin;
    k = 6;
    goto L10083;
L10081:
    --k;
L10083:
    if (-(k - 1) > 0) {
	goto L10082;
    }
    if (cvrss[k - 1] > cvmin) {
	goto L10101;
    }
    goto L10082;
L10101:
    goto L10081;
L10082:
    *span = cvspan[k - 1];
L10071:
L10021:
    smth_(&x[1], &y[1], &w[1], span, dof, n, &c__0, &smo[1], &s0, rss, &scrat[
	    1]);
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	smo[i] += s0;
/* L10111: */
    }
/* L10112: */
    return 0;
} /* rlsmo_ */

/* Subroutine */ int smth_(x, y, w, span, dof, n, cross, smo, s0, rss, scrat)
doublereal *x, *y, *w, *span;
real *dof;
integer *n, *cross;
doublereal *smo;
real *s0;
doublereal *rss, *scrat;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer line;
    static doublereal xbar, ybar;
    static integer ntie;
    static doublereal sumw;
    static real xout, yout;
    static integer i, j;
    static real r;
    static integer ibold, ibnew, ispan, itold, itnew, m0, jj;
    static real wt;
    static integer fixeds, is2;
    static doublereal cov, var;
    static real win, xin, yin;

    /* Parameter adjustments */
    --scrat;
    --smo;
    --w;
    --y;
    --x;

    /* Function Body */
    line = 1;
    fixeds = 1;
    if (*span >= (float)1.) {
	goto L10131;
    }
    line = 0;
L10131:
    xbar = x[1];
    ybar = y[1];
    cov = (float)0.;
    var = (float)0.;
    sumw = w[1];
    if (line != 1) {
	goto L10151;
    }
    i__1 = *n;
    for (i = 2; i <= i__1; ++i) {
	xin = x[i];
	yin = y[i];
	win = w[i];
	xbar = (sumw * xbar + xin * win) / (sumw + win);
	ybar = (sumw * ybar + yin * win) / (sumw + win);
	cov += win * (xin - xbar) * (yin - ybar) * (sumw + win) / sumw;
/* Computing 2nd power */
	d__1 = xin - xbar;
	var += win * (d__1 * d__1) * (sumw + win) / sumw;
	sumw += win;
/* L10161: */
    }
/* L10162: */
    i = 1;
    goto L10173;
L10171:
    ++i;
L10173:
    if (i > *n) {
	goto L10172;
    }
    if (! (*cross == 1)) {
	goto L10191;
    }
    xout = x[i];
    yout = y[i];
    win = w[i];
    cov -= win * (xout - xbar) * (yout - ybar) * sumw / (sumw - win);
/* Computing 2nd power */
    d__1 = xout - xbar;
    var -= win * (d__1 * d__1) * sumw / (sumw - win);
    xbar = (sumw * xbar - win * xout) / (sumw - win);
    ybar = (sumw * ybar - win * yout) / (sumw - win);
    sumw -= win;
L10191:
    if (var <= (float)0.) {
	goto L10211;
    }
    smo[i] = cov * (x[i] - xbar) / var;
    goto L10221;
L10211:
    smo[i] = 0.;
L10221:
/* L10201: */
    if (! (*cross == 1)) {
	goto L10241;
    }
    xin = x[i];
    yin = y[i];
    win = w[i];
    xbar = (sumw * xbar + xin * win) / (sumw + win);
    ybar = (sumw * ybar + yin * win) / (sumw + win);
    cov += win * (xin - xbar) * (yin - ybar) * (sumw + win) / sumw;
/* Computing 2nd power */
    d__1 = xin - xbar;
    var += win * (d__1 * d__1) * (sumw + win) / sumw;
    sumw += win;
L10241:
    goto L10171;
L10172:
    *s0 = ybar;
    scrat[1] = cov / var;
    *dof = (float)1.;
    goto L10251;
L10151:
    itold = 1;
    ibold = 1;
    *dof = (float)-1.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	scrat[i] = y[i];
/* L10261: */
    }
/* L10262: */
    if (! (*cross == 0)) {
	goto L10281;
    }
    i = 0;
L10291:
    if (i >= *n - 1) {
	goto L10292;
    }
    ++i;
    m0 = i;
L10301:
    if (x[i + 1] > x[i]) {
	goto L10302;
    }
    ++i;
    if (i < *n) {
	goto L10301;
    }
L10302:
    if (i == m0) {
	goto L10291;
    }
    ntie = i - m0 + 1;
    r = (float)0.;
    wt = (float)0.;
    i__1 = i;
    for (jj = m0; jj <= i__1; ++jj) {
	j = jj;
	r += y[j] * w[j];
	wt += w[j];
/* L10311: */
    }
/* L10312: */
    r /= wt;
    i__1 = i;
    for (j = m0; j <= i__1; ++j) {
	y[j] = r;
/* L10321: */
    }
/* L10322: */
    goto L10291;
L10292:
L10281:
    ispan = (integer) (*n * *span);
    if (! (fixeds == 1)) {
	goto L10341;
    }
    is2 = ispan / 2;
    if (is2 >= 1) {
	goto L10361;
    }
    is2 = 1;
L10361:
L10341:
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
/* Computing MIN */
	i__2 = i + is2;
	itnew = min(i__2,*n);
/* Computing MAX */
	i__2 = i - is2;
	ibnew = max(i__2,1);
L10381:
	if (itold >= itnew) {
	    goto L10382;
	}
	++itold;
	xin = x[itold];
	yin = y[itold];
	win = w[itold];
	xbar = (sumw * xbar + xin * win) / (sumw + win);
	ybar = (sumw * ybar + yin * win) / (sumw + win);
	cov += win * (xin - xbar) * (yin - ybar) * (sumw + win) / sumw;
/* Computing 2nd power */
	d__1 = xin - xbar;
	var += win * (d__1 * d__1) * (sumw + win) / sumw;
	sumw += win;
	goto L10381;
L10382:
L10391:
	if (ibold <= ibnew) {
	    goto L10392;
	}
	--ibold;
	xin = x[ibold];
	yin = y[ibold];
	win = w[ibold];
	xbar = (sumw * xbar + xin * win) / (sumw + win);
	ybar = (sumw * ybar + yin * win) / (sumw + win);
	cov += win * (xin - xbar) * (yin - ybar) * (sumw + win) / sumw;
/* Computing 2nd power */
	d__1 = xin - xbar;
	var += win * (d__1 * d__1) * (sumw + win) / sumw;
	sumw += win;
	goto L10391;
L10392:
L10401:
	if (itold <= itnew) {
	    goto L10402;
	}
	xout = x[itold];
	yout = y[itold];
	win = w[itold];
	cov -= win * (xout - xbar) * (yout - ybar) * sumw / (sumw - win);
/* Computing 2nd power */
	d__1 = xout - xbar;
	var -= win * (d__1 * d__1) * sumw / (sumw - win);
	xbar = (sumw * xbar - win * xout) / (sumw - win);
	ybar = (sumw * ybar - win * yout) / (sumw - win);
	sumw -= win;
	--itold;
	goto L10401;
L10402:
L10411:
	if (ibold >= ibnew) {
	    goto L10412;
	}
	xout = x[ibold];
	yout = y[ibold];
	win = w[ibold];
	cov -= win * (xout - xbar) * (yout - ybar) * sumw / (sumw - win);
/* Computing 2nd power */
	d__1 = xout - xbar;
	var -= win * (d__1 * d__1) * sumw / (sumw - win);
	xbar = (sumw * xbar - win * xout) / (sumw - win);
	ybar = (sumw * ybar - win * yout) / (sumw - win);
	sumw -= win;
	++ibold;
	goto L10411;
L10412:
	if (! (*cross == 1)) {
	    goto L10431;
	}
	xout = x[i];
	yout = y[i];
	win = w[i];
	cov -= win * (xout - xbar) * (yout - ybar) * sumw / (sumw - win);
/* Computing 2nd power */
	d__1 = xout - xbar;
	var -= win * (d__1 * d__1) * sumw / (sumw - win);
	xbar = (sumw * xbar - win * xout) / (sumw - win);
	ybar = (sumw * ybar - win * yout) / (sumw - win);
	sumw -= win;
L10431:
	if (var <= (float)0.) {
	    goto L10451;
	}
	smo[i] = ybar + cov * (x[i] - xbar) / var;
/* Computing 2nd power */
	d__1 = x[i] - xbar;
	*dof = *dof + w[i] / sumw + w[i] * (d__1 * d__1) / var;
	goto L10461;
L10451:
	smo[i] = ybar;
	*dof += w[i] / sumw;
L10461:
/* L10441: */
	if (! (*cross == 1)) {
	    goto L10481;
	}
	xin = x[i];
	yin = y[i];
	win = w[i];
	xbar = (sumw * xbar + xin * win) / (sumw + win);
	ybar = (sumw * ybar + yin * win) / (sumw + win);
	cov += win * (xin - xbar) * (yin - ybar) * (sumw + win) / sumw;
/* Computing 2nd power */
	d__1 = xin - xbar;
	var += win * (d__1 * d__1) * (sumw + win) / sumw;
	sumw += win;
L10481:
/* L10371: */
	;
    }
/* L10372: */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	y[i] = scrat[i];
/* L10491: */
    }
/* L10492: */
    if (*cross != 0) {
	goto L10511;
    }
    i = 0;
L10521:
    if (i >= *n - 1) {
	goto L10522;
    }
    ++i;
    m0 = i;
L10531:
    if (x[i + 1] > x[i]) {
	goto L10532;
    }
    ++i;
    if (i < *n) {
	goto L10531;
    }
L10532:
    if (i == m0) {
	goto L10521;
    }
    ntie = i - m0 + 1;
    r = (float)0.;
    wt = (float)0.;
    i__1 = i;
    for (jj = m0; jj <= i__1; ++jj) {
	j = jj;
	r += smo[j] * w[j];
	wt += w[j];
/* L10541: */
    }
/* L10542: */
    r /= wt;
    i__1 = i;
    for (j = m0; j <= i__1; ++j) {
	smo[j] = r;
/* L10551: */
    }
/* L10552: */
    goto L10521;
L10522:
L10511:
    ybar = (float)0.;
    sumw = (float)0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	ybar += w[i] * y[i];
	sumw += w[i];
/* L10561: */
    }
/* L10562: */
    ybar /= sumw;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	smo[i] -= ybar;
/* L10571: */
    }
/* L10572: */
    *s0 = ybar;
L10251:
/* L10141: */
    *rss = (float)0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
/* Computing 2nd power */
	d__1 = y[i] - *s0 - smo[i];
	*rss += w[i] / sumw * (d__1 * d__1);
/* L10581: */
    }
/* L10582: */
    return 0;
} /* smth_ */

