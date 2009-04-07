/* Kenneth Clarkson's hull routines merged in a single file: 
   by R. Wenger, November, 2001 */
/* clarkson_convex_hull: interface routine for simple retrieval 
   of convex hull simplices (and their vertices)
   by R. Wenger, November, 2001 */

/*
 * Ken Clarkson wrote this.  Copyright (c) 1995 by AT&T..
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */

/*
  Includes Clarkson's modification of:

           An implementation of top-down splaying with sizes
             D. Sleator <sleator@cs.cmu.edu>, January 1994.
*/


#include "clarkson_hull.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* hull.c */

site p;
long pnum;

int	rdim,	/* region dimension: (max) number of sites specifying region */
	cdim,	/* number of sites currently specifying region */
	site_size, /* size of malloc needed for a site */
	point_size;  /* size of malloc needed for a point */

STORAGE(simplex)

#define push(x) *(st+tms++) = x;
#define pop(x)  x = *(st + --tms);

void *visit_triang_gen(simplex *s, visit_func *visit, test_func *test) {
/* 
 * starting at s, visit simplices t such that test(s,i,0) is true,
 * and t is the i'th neighbor of s;
 * apply visit function to all visited simplices;
 * when visit returns nonNULL, exit and return its value
 */
	neighbor *sn;
	void *v;
	simplex *t;
	int i;
	long tms = 0;
	static long vnum = -1;
	static long ss = 2000;
	static simplex **st;

	vnum--;
	if (!st) { 
	  /* bug fix: 10-17-2006 by R. Wenger */
	  st=(simplex**)malloc((ss+MAXDIM+1)*sizeof(simplex*));
	  assert(st != 0);
	};

	if (s) push(s);

	while (tms) {

		if (tms>ss) {
		  DEBEXP(-1,tms);
		  /* bug fix: 10-17-2006 by R. Wenger */
		  st=(simplex**)  
		    realloc(st, ((ss+=ss)+MAXDIM+1)*sizeof(simplex*));
		  assert(st != 0);
		};
		pop(t);
		if (!t || t->visit == vnum) continue;
		t->visit = vnum;
		if (v=(*visit)(t,0)) {return v;}
		for (i=-1,sn = t->neigh-1;i<cdim;i++,sn++)
			if ((sn->simp->visit != vnum) && sn->simp && test(t,i,0))
				push(sn->simp);
	}

	return NULL;
}

int truet(simplex *s, int i, void *dum) {return 1;}

void *visit_triang(simplex *root, visit_func *visit)
/* visit the whole triangulation */
	{return visit_triang_gen(root, visit, truet);}


int hullt(simplex *s, int i, void *dummy) {return i>-1;}

void *facet_test(simplex *s, void *dummy) {return (!s->peak.vert) ? s : NULL;}

void *visit_hull(simplex *root, visit_func *visit)
/* visit all simplices with facets of the current hull */
	{return visit_triang_gen(visit_triang(root, &facet_test),
			visit, hullt);}



#define lookup(a,b,what,whatt)						\
{									\
	int i;								\
	neighbor *x;							\
	for (i=0, x = a->neigh; (x->what != b) && (i<cdim) ; i++, x++);	\
	if (i<cdim)							\
		return x;						\
	else {								\
		fprintf(DFILE,"adjacency failure,op_" #what ":\n");	\
		DEBTR(-10)						\
		print_simplex_f(a, DFILE, &print_neighbor_full);	\
		print_##whatt(b, DFILE);				\
		fprintf(DFILE,"---------------------\n");		\
		print_triang(a,DFILE, &print_neighbor_full);		\
		exit(1);						\
		return 0;						\
	}								\
}									\


neighbor *op_simp(simplex *a, simplex *b) {lookup(a,b,simp,simplex)}
	/* the neighbor entry of a containing b */

neighbor *op_vert(simplex *a, site b) {lookup(a,b,vert,site)}
	/* the neighbor entry of a containing b */


void connect(simplex *s) {
/* make neighbor connections between newly created simplices incident to p */

	site xf,xb,xfi;
	simplex *sb, *sf, *seen;
	int i;
	neighbor *sn;

	if (!s) return;
	assert(!s->peak.vert
		&& s->peak.simp->peak.vert==p
		&& !op_vert(s,p)->simp->peak.vert);
	if (s->visit==pnum) return;
	s->visit = pnum;
	seen = s->peak.simp;
	xfi = op_simp(seen,s)->vert;
	for (i=0, sn = s->neigh; i<cdim; i++,sn++) {
		xb = sn->vert;
		if (p == xb) continue;
		sb = seen;
		sf = sn->simp;
		xf = xfi;
		if (!sf->peak.vert) {	/* are we done already? */
			sf = op_vert(seen,xb)->simp;
			if (sf->peak.vert) continue;				
		} else do {
			xb = xf;
			xf = op_simp(sf,sb)->vert;
			sb = sf;
			sf = op_vert(sb,xb)->simp;
		} while (sf->peak.vert);

		sn->simp = sf;
		op_vert(sf,xf)->simp = s;

		connect(sf);
	}

}


				
simplex *make_facets(simplex *seen) {
/*
 * visit simplices s with sees(p,s), and make a facet for every neighbor
 * of s not seen by p
 */

	simplex *n;
	static simplex *ns;
	neighbor *bn;
	int i;


	if (!seen) return NULL;
	DEBS(-1) assert(sees(p,seen) && !seen->peak.vert); EDEBS
	seen->peak.vert = p;

	for (i=0,bn = seen->neigh; i<cdim; i++,bn++) {
		n = bn->simp;
		if (pnum != n->visit) {
			n->visit = pnum;
			if (sees(p,n)) make_facets(n);
		} 
		if (n->peak.vert) continue;
		copy_simp(ns,seen);
		ns->visit = 0;
		ns->peak.vert = 0;
		ns->normal = 0;
		ns->peak.simp = seen;
/*		ns->Sb -= ns->neigh[i].basis->sqb; */
		NULLIFY(basis_s,ns->neigh[i].basis);
		ns->neigh[i].vert = p;
		bn->simp = op_simp(n,seen)->simp = ns;
	}
	return ns;
}



simplex *extend_simplices(simplex *s) {
/*
 * p lies outside flat containing previous sites;
 * make p a vertex of every current simplex, and create some new simplices
 */

	int	i,
		ocdim=cdim-1;
	simplex *ns;
	neighbor *nsn;

	if (s->visit == pnum) return s->peak.vert ? s->neigh[ocdim].simp : s;
	s->visit = pnum;
	s->neigh[ocdim].vert = p;
	NULLIFY(basis_s,s->normal);
	NULLIFY(basis_s,s->neigh[0].basis);
	if (!s->peak.vert) {
		s->neigh[ocdim].simp = extend_simplices(s->peak.simp);
		return s;
	} else {
		copy_simp(ns,s);
		s->neigh[ocdim].simp = ns;
		ns->peak.vert = NULL;
		ns->peak.simp = s;
		ns->neigh[ocdim] = s->peak;
		inc_ref(basis_s,s->peak.basis);
		for (i=0,nsn=ns->neigh;i<cdim;i++,nsn++)
			nsn->simp = extend_simplices(nsn->simp);
	}
	return ns;
}


simplex *search(simplex *root) {
/* return a simplex s that corresponds to a facet of the 
 * current hull, and sees(p, s) */

	simplex *s;
	static simplex **st;
	static long ss = MAXDIM;
	neighbor *sn;
	int i;
	long tms = 0;

	if (!st) st = (simplex **)malloc((ss+MAXDIM+1)*sizeof(simplex*));
	push(root->peak.simp);
	root->visit = pnum;
	if (!sees(p,root))
		for (i=0,sn=root->neigh;i<cdim;i++,sn++) push(sn->simp);
	while (tms) {
	  if (tms>ss) {
	    /* bug fix: 10-17-2006 by R. Wenger */
	    st=(simplex**)
	      realloc(st, ((ss+=ss)+MAXDIM+1)*sizeof(simplex*));
	    assert(st != 0);
	  };

	  pop(s);
	  if (s->visit == pnum) continue;
	  s->visit = pnum;
	  if (!sees(p,s)) continue;
	  if (!s->peak.vert) return s;
	  for (i=0, sn=s->neigh; i<cdim; i++,sn++) push(sn->simp);
	}
	return NULL;
}


point get_another_site(void) {

/*	static int scount =0; */
	point pnext;

/*	if (!(++scount%1000)) {fprintf(DFILE,"site %d...", scount);} */
/*	check_triang(); */
	pnext = (*get_site)();
	if (!pnext) return NULL;
	pnum = site_num(pnext)+2;
	return pnext;
}



void buildhull (simplex *root) {


	while (cdim < rdim) {
		p = get_another_site();
		if (!p) return;
		if (out_of_flat(root,p)) {
		  extend_simplices(root);
		}
		else {
		  connect(make_facets(search(root)));
		}
	}

	while (p = get_another_site()) {
	  connect(make_facets(search(root)));
	}
}

/* ch.c */

short check_overshoot_f=0;


simplex *ch_root;

#define NEARZERO(d)	((d) < FLT_EPSILON && (d) > -FLT_EPSILON)
#define SMALL (100*FLT_EPSILON)*(100*FLT_EPSILON)

#define SWAP(X,a,b) {X t; t = a; a = b; b = t;}

#define DMAX 

double Huge;

#define check_overshoot(x)							\
	{if (CHECK_OVERSHOOT && check_overshoot_f && ((x)>9e15))		\
		warning(-20, overshot exact arithmetic)}			\


#define DELIFT 0
int basis_vec_size;


#define lookupshort(a,b,whatb,c,whatc)					\
{									\
	int i;								\
	neighbor *x;							\
	c = NULL;							\
	for (i=-1, x = a->neigh-1; (x->whatb != b) && (i<cdim) ; i++, x++);\
	if (i<cdim) c = x->whatc;					\
}									\

	
Coord Vec_dot(point x, point y) {
	int i;
	Coord sum = 0;
	for (i=0;i<rdim;i++) sum += x[i] * y[i];
	return sum;
}

Coord Vec_dot_pdim(point x, point y) {
	int i;
	Coord sum = 0;
	for (i=0;i<pdim;i++) sum += x[i] * y[i];
/*	check_overshoot(sum); */
	return sum;
}

Coord Norm2(point x) {
	int i;
	Coord sum = 0;
	for (i=0;i<rdim;i++) sum += x[i] * x[i];
	return sum;
}

void Ax_plus_y(Coord a, point x, point y) {
	int i;
	for (i=0;i<rdim;i++) {
		*y++ += a * *x++;
	}
}

void Ax_plus_y_test(Coord a, point x, point y) {
	int i;
	for (i=0;i<rdim;i++) {
		check_overshoot(*y + a * *x);
		*y++ += a * *x++;
	}
}

void Vec_scale(int n, Coord a, Coord *x)
{
    register Coord *xx = x,
		*xend = xx + n;
	while (xx!=xend) *xx++ *= a;
}

void Vec_scale_test(int n, Coord a, Coord *x)
{
    register Coord *xx = x,
		*xend = xx + n	;

	check_overshoot(a);

	while (xx!=xend) {
		*xx *= a;
		check_overshoot(*xx);
		xx++;
	}
}



int exact_bits;
float b_err_min, b_err_min_sq;

double logb(double); /* on SGI machines: returns floor of log base 2 */

static short vd;
static basis_s	tt_basis = {0,1,-1,0,0,0},
		*tt_basisp = &tt_basis,
		*infinity_basis;


STORAGE(basis_s)

typedef Coord site_struct;

Coord	infinity[10]={57.2,0,0,0,0}; /* point at infinity for vd; value not used */
	
void print_site(site p, FILE* F)
	{print_point(F, pdim,p);fprintf(F, "\n");}

#define VA(x) ((x)->vecs+rdim)
#define VB(x) ((x)->vecs)




/* tables for runtime stats */
int A[100]={0}, B[100] ={0}, C[100] = {0}, D[100] ={0};
int tot =0, totinf=0, bigt=0; 

#define two_to(x) ( ((x)<20) ? 1<<(x) : ldexp(1,(x)) )



double sc(basis_s *v,simplex *s, int k, int j) {
/* amount by which to scale up vector, for reduce_inner */

	double		labound;
	static int	lscale;
	static double	max_scale,
			ldetbound,
			Sb;

	if (j<10) {
		labound = logb(v->sqa)/2;
		max_scale = exact_bits - labound - 0.66*(k-2)-1  -DELIFT;
		if (max_scale<1) {
			warning(-10, overshot exact arithmetic);
			max_scale = 1;

		}

		if (j==0) {
			int	i;
			neighbor *sni;
			basis_s *snib;

			ldetbound = DELIFT;

			Sb = 0;
			for (i=k-1,sni=s->neigh+k-1;i>0;i--,sni--) {
				snib = sni->basis;
				Sb += snib->sqb;
				ldetbound += logb(snib->sqb)/2 + 1;
				ldetbound -= snib->lscale;
			}
		}
	}
	if (ldetbound - v->lscale + logb(v->sqb)/2 + 1 < 0) {
		DEBS(-2)
			DEBTR(-2) DEBEXP(-2, ldetbound)
			print_simplex_f(s, DFILE, &print_neighbor_full);
			print_basis(DFILE,v);
		EDEBS			
		return 0;					
	} else {
		lscale = logb(2*Sb/(v->sqb + v->sqa*b_err_min))/2;	
		if (lscale > max_scale) {
			lscale = max_scale;
		} else if (lscale<0) lscale = 0;
		v->lscale += lscale;
		return two_to(lscale);
	}
}


double lower_terms(basis_s* v) {

	point vp = v->vecs;
	int i,j,h,hh=0;
	int facs[6] = {2,3,5,7,11,13};
	double out = 1;

DEBTR(-10) print_basis(DFILE, v); printf("\n");
DEBTR(0)

	for (j=0;j<6;j++) do {
		for (i=0; i<2*rdim && facs[j]*floor(vp[i]/facs[j])==vp[i];i++);
		if (h = (i==2*rdim)) {
			hh=1;
			out *= facs[j];
			for (i=0;i<2*rdim; i++) vp[i]/=facs[j];
		}
	} while (h);
/*	if (hh) {DEBTR(-10)  print_basis(DFILE, v);} */
	return out;
}

double lower_terms_point(point vp) {

	int i,j,h,hh=0;
	int facs[6] = {2,3,5,7,11,13};
	double out = 1;

	for (j=0;j<6;j++) do {
		for (i=0; i<2*rdim && facs[j]*floor(vp[i]/facs[j])==vp[i];i++);
		if (h = (i==2*rdim)) {
			hh=1;
			out *= facs[j];
			for (i=0;i<2*rdim; i++) vp[i]/=facs[j];
		}
	} while (h);
	return out;
}


int reduce_inner(basis_s *v, simplex *s, int k) {

	point	va = VA(v),
		vb = VB(v);
	int	i,j;
	double	dd,
		ldetbound = 0,
		Sb = 0;
	double	scale;
	basis_s	*snibv;
	neighbor *sni;
	static int failcount;

/*	lower_terms(v); */
	v->sqa = v->sqb = Norm2(vb);
	if (k<=1) {
		memcpy(vb,va,basis_vec_size);
		return 1;
	}
/*	if (vd) {
		snibv = s->neigh[1].basis;
		scale = floor(sqrt(snibv->sqa/v->sqa));
		if (scale > 1) Vec_scale(rdim,scale,va);
		dd = Vec_dot(VA(snibv),va)/snibv->sqa;
		dd = -floor(0.5+dd);
		Ax_plus_y( dd, VA(snibv), va);
	}
*/		
	for (j=0;j<250;j++) {

		memcpy(vb,va,basis_vec_size);
		for (i=k-1,sni=s->neigh+k-1;i>0;i--,sni--) {
			snibv = sni->basis;
			dd = -Vec_dot(VB(snibv),vb)/ snibv->sqb;
			Ax_plus_y( dd, VA(snibv), vb);		
		}
		v->sqb = Norm2(vb);		
		v->sqa = Norm2(va);
		
		if (2*v->sqb >= v->sqa) {B[j]++; return 1;}

		Vec_scale_test(rdim,scale = sc(v,s,k,j),va);
		
		for (i=k-1,sni=s->neigh+k-1;i>0;i--,sni--) {
			snibv = sni->basis;
			dd = -Vec_dot(VB(snibv),va)/snibv->sqb;	
			dd = floor(dd+0.5);
			Ax_plus_y_test( dd, VA(snibv), va);
		}		
	}
	if (failcount++<10) {
		DEB(-8, reduce_inner failed on:)
		DEBTR(-8) print_basis(DFILE, v); 
		print_simplex_f(s, DFILE, &print_neighbor_full);
	}
	return 0;
}

#define trans(z,p,q) {int i; for (i=0;i<pdim;i++) z[i+rdim] = z[i] = p[i] - q[i];}
#define lift(z,s) {if (vd) z[2*rdim-1] =z[rdim-1]= ldexp(Vec_dot_pdim(z,z), -DELIFT);}
				/*not scaling lift to 2^-DELIFT */



int reduce(basis_s **v, point p, simplex *s, int k) {

	point	z;
	point	tt = s->neigh[0].vert;

	if (!*v) NEWLRC(basis_s,(*v))
	else (*v)->lscale = 0;
	z = VB(*v);
	if (vd) {
		if (p==infinity) memcpy(*v,infinity_basis,basis_s_size);
		else {trans(z,p,tt); lift(z,s);}
	} else trans(z,p,tt);
	return reduce_inner(*v,s,k);
}


void get_basis_sede(simplex *s) {

	int	k=1;
	neighbor *sn = s->neigh+1,
		 *sn0 = s->neigh;

	if (vd && sn0->vert == infinity && cdim >1) {
		SWAP(neighbor, *sn0, *sn );
		NULLIFY(basis_s,sn0->basis);
		sn0->basis = tt_basisp;
		tt_basisp->ref_count++;
	} else {
		if (!sn0->basis) {
			sn0->basis = tt_basisp;
			tt_basisp->ref_count++;
		} else while (k < cdim && sn->basis) {k++;sn++;}
	}
	while (k<cdim) {
		NULLIFY(basis_s,sn->basis);
		reduce(&sn->basis,sn->vert,s,k);
		k++; sn++;
	}
}


int out_of_flat(simplex *root, point p) {

	static neighbor p_neigh={0,0,0};

	if (!p_neigh.basis) p_neigh.basis = (basis_s*) malloc(basis_s_size);

	p_neigh.vert = p;
	cdim++;
	root->neigh[cdim-1].vert = root->peak.vert;
	NULLIFY(basis_s,root->neigh[cdim-1].basis);
	get_basis_sede(root);
	if (vd && root->neigh[0].vert == infinity) return 1;
	reduce(&p_neigh.basis,p,root,cdim);
	if (p_neigh.basis->sqa != 0) return 1;
	cdim--;
	return 0;
}


double cosangle_sq(basis_s* v,basis_s* w)  {
	double dd;
	point	vv=v->vecs,
		wv=w->vecs;
	dd = Vec_dot(vv,wv);
	return dd*dd/Norm2(vv)/Norm2(wv);
}


int check_perps(simplex *s) {

	static basis_s *b = NULL;
	point 	z,y;
	point	tt;
	double	dd;
	int i,j;

	for (i=1; i<cdim; i++) if (NEARZERO(s->neigh[i].basis->sqb)) return 0;
	if (!b) {
	  /* bug fix: 10-17-2006 by R. Wenger */
	  b = (basis_s*)malloc(basis_s_size); 
	  assert(b != 0);
	};
	b->lscale = 0;
	z = VB(b);
	tt = s->neigh[0].vert;
	for (i=1;i<cdim;i++) {
		y = s->neigh[i].vert;
		if (vd && y==infinity) memcpy(b, infinity_basis, basis_s_size);
		else {trans(z,y,tt); lift(z,s);}
		if (s->normal && cosangle_sq(b,s->normal)>b_err_min_sq) {DEBS(0)
			DEB(0,bad normal) DEBEXP(0,i) DEBEXP(0,dd)
			print_simplex_f(s, DFILE, &print_neighbor_full);
			EDEBS
			return 0;
		}
		for (j=i+1;j<cdim;j++) {
			if (cosangle_sq(b,s->neigh[j].basis)>b_err_min_sq) {
				DEBS(0)
				DEB(0,bad basis)DEBEXP(0,i) DEBEXP(0,j) DEBEXP(0,dd)
				DEBTR(-8) print_basis(DFILE, b);
				print_simplex_f(s, DFILE, &print_neighbor_full);
				EDEBS
				return 0;
			}
		}
	}
	return 1;
}


void get_normal_sede(simplex *s) {

	neighbor *rn;
	int i,j;

	get_basis_sede(s);
	if (rdim==3 && cdim==3) {
		point	c,
			a = VB(s->neigh[1].basis),
			b = VB(s->neigh[2].basis);
		NEWLRC(basis_s,s->normal);
		c = VB(s->normal);
		c[0] = a[1]*b[2] - a[2]*b[1];
		c[1] = a[2]*b[0] - a[0]*b[2];
		c[2] = a[0]*b[1] - a[1]*b[0];
		s->normal->sqb = Norm2(c);
		for (i=cdim+1,rn = ch_root->neigh+cdim-1; i; i--, rn--) {
			for (j = 0; j<cdim && rn->vert != s->neigh[j].vert;j++);
			if (j<cdim) continue;
			if (rn->vert==infinity) {
				if (c[2] > -b_err_min) continue;
			} else  if (!sees(rn->vert,s)) continue;
			c[0] = -c[0]; c[1] = -c[1]; c[2] = -c[2];
			break;
		}
		DEBS(-1) if (!check_perps(s)) exit(1); EDEBS
		return;
	}	
		
	for (i=cdim+1,rn = ch_root->neigh+cdim-1; i; i--, rn--) {
		for (j = 0; j<cdim && rn->vert != s->neigh[j].vert;j++);
		if (j<cdim) continue;
		reduce(&s->normal,rn->vert,s,cdim);
		if (s->normal->sqb != 0) break;
	}

	DEBS(-1) if (!check_perps(s)) {DEBTR(-1) exit(1);} EDEBS

}


void get_normal(simplex *s) {get_normal_sede(s); return;}

int sees(site p, simplex *s) {

	static basis_s *b = NULL;
	point	tt,zz;
	double	dd,dds;
	int i;


	if (!b) {
	  /* bug fix: 10-17-2006 by R. Wenger */
	  b = (basis_s*)malloc(basis_s_size);
	  assert(b != 0);
	};
	b->lscale = 0;
	zz = VB(b);
	if (cdim==0) return 0;
	if (!s->normal) {
		get_normal_sede(s);
		for (i=0;i<cdim;i++) NULLIFY(basis_s,s->neigh[i].basis);
	}
	tt = s->neigh[0].vert;
	if (vd) {
		if (p==infinity) memcpy(b,infinity_basis,basis_s_size);
		else {trans(zz,p,tt); lift(zz,s);}
	} else trans(zz,p,tt);
	for (i=0;i<3;i++) {
		dd = Vec_dot(zz,s->normal->vecs);	
		if (dd == 0.0) {
			DEBS(-7) DEB(-6,degeneracy:); DEBEXP(-6,site_num(p));
			print_site(p, DFILE); print_simplex_f(s, DFILE, &print_neighbor_full); EDEBS
			return 0;
		} 
		dds = dd*dd/s->normal->sqb/Norm2(zz);
		if (dds > b_err_min_sq) return (dd<0);
		get_basis_sede(s);
		reduce_inner(b,s,cdim);
	}
	DEBS(-7) if (i==3) {
		DEB(-6, looped too much in sees);
		DEBEXP(-6,dd) DEBEXP(-6,dds) DEBEXP(-6,site_num(p));
		print_simplex_f(s, DFILE, &print_neighbor_full); exit(1);}
	EDEBS
	return 0;
}





double radsq(simplex *s) {

	Coord z[MAXDIM];
	point n;
	neighbor *sn;
	int i;

		/* square of ratio of circumcircle radius to
			max edge length for Delaunay tetrahedra */

	
	for (i=0,sn=s->neigh;i<cdim;i++,sn++)
		if (sn->vert == infinity) return Huge;

	if (!s->normal) get_normal_sede(s);

					/* compute circumradius */
	n = s->normal->vecs;

	if (NEARZERO(n[rdim-1])) {return Huge;}

	return Vec_dot_pdim(n,n)/4/n[rdim-1]/n[rdim-1];
}


void *zero_marks(simplex * s, void *dum) { s->mark = 0; return NULL; }

void *one_marks(simplex * s, void *dum) {s->mark = 1; return NULL;}

void *show_marks(simplex * s, void *dum) {printf("%d",s->mark); return NULL;}


#define swap_points(a,b) {point t; t=a; a=b; b=t;}

int alph_test(simplex *s, int i, void *alphap) {
				/*returns 1 if not an alpha-facet */

	simplex *si;
	double	rs,rsi,rsfi;
	neighbor *scn, *sin;
	int k, nsees, ssees;
	static double alpha;

	if (alphap) {alpha=*(double*)alphap; if (!s) return 1;}
	if (i==-1) return 0;

	si = s->neigh[i].simp;
	scn = s->neigh+cdim-1;
	sin = s->neigh+i;
	nsees = 0;

	for (k=0;k<cdim;k++) if (s->neigh[k].vert==infinity && k!=i) return 1;
	rs = radsq(s);
	rsi = radsq(si);

	if (rs < alpha &&  rsi < alpha) return 1;

	swap_points(scn->vert,sin->vert);
	NULLIFY(basis_s, s->neigh[i].basis);
	cdim--;
	get_basis_sede(s);
	reduce(&s->normal,infinity,s,cdim);
	rsfi = radsq(s);

	for (k=0;k<cdim;k++) if (si->neigh[k].simp==s) break;

	ssees = sees(scn->vert,s);
	if (!ssees) nsees = sees(si->neigh[k].vert,s);
	swap_points(scn->vert,sin->vert);
	cdim++;
	NULLIFY(basis_s, s->normal);
	NULLIFY(basis_s, s->neigh[i].basis);

	if (ssees) return alpha<rs;
	if (nsees) return alpha<rsi;

	assert(rsfi<=rs+FLT_EPSILON && rsfi<=rsi+FLT_EPSILON);

	return alpha<=rsfi;
}


void *conv_facetv(simplex *s, void *dum) {
	int i;
	for (i=0;i<cdim;i++) if (s->neigh[i].vert==infinity) {return s;}
	return NULL;
}

short mi[MAXPOINTS], mo[MAXPOINTS];

void *mark_points(simplex *s, void *dum) {
	int i, snum;
	neighbor *sn;

	for  (i=0,sn=s->neigh;i<cdim;i++,sn++) {
		if (sn->vert==infinity) continue;
		snum = site_num(sn->vert);
		if (s->mark) mo[snum] = 1;
		else mi[snum] = 1;
	}
	return NULL;
}

void* visit_outside_ashape(simplex *root, visit_func visit) {
	return visit_triang_gen(visit_hull(root, conv_facetv), visit, alph_test);
}

int check_ashape(simplex *root, double alpha) {

	int i;

	for (i=0;i<MAXPOINTS;i++) {mi[i] = mo[i] = 0;}

	visit_hull(root, zero_marks);

	alph_test(0,0,&alpha);
	visit_outside_ashape(root, one_marks);

	visit_hull(root, mark_points);

	for (i=0;i<MAXPOINTS;i++) if (mo[i] && !mi[i]) {return 0;}

	return 1;
}

double find_alpha(simplex *root) {

	int i;
	float al=0,ah,am;

	for (ah=i=0;i<pdim;i++) ah += (maxs[i]-mins[i])*(maxs[i]-mins[i]);
	assert(check_ashape(root,ah));
	for (i=0;i<17;i++) {
		if (check_ashape(root, am = (al+ah)/2)) ah = am;
		else al = am;
		if ((ah-al)/ah<.5) break;
	}
	return 1.1*ah;
}






void vols(fg *f, Tree *t, basis_s* n, int depth) {

	static simplex *s;
	static neighbor *sn;
	int tdim = cdim;
	basis_s *nn = 0;
	int signum;
	point nnv;
	double sqq;


	if (!t) {return;}

	if (!s) {NEWL(simplex,s); sn = s->neigh;}
	cdim = depth;
	s->normal = n;
	if (depth>1 && sees(t->key,s)) signum = -1; else signum = 1;
	cdim = tdim;

	if (t->fgs->dist == 0) {
		sn[depth-1].vert = t->key;
		NULLIFY(basis_s,sn[depth-1].basis);
		cdim = depth; get_basis_sede(s); cdim = tdim;
		reduce(&nn, infinity, s, depth);
		nnv = nn->vecs;
		if (t->key==infinity || f->dist==Huge || NEARZERO(nnv[rdim-1]))
			t->fgs->dist = Huge;
		else
			t->fgs->dist = Vec_dot_pdim(nnv,nnv)
						/4/nnv[rdim-1]/nnv[rdim-1];
		if (!t->fgs->facets) t->fgs->vol = 1;
		else vols(t->fgs, t->fgs->facets, nn, depth+1);
	}

	assert(f->dist!=Huge || t->fgs->dist==Huge);
	if (t->fgs->dist==Huge || t->fgs->vol==Huge) f->vol = Huge;
	else {
		sqq = t->fgs->dist - f->dist;
		if (NEARZERO(sqq)) f->vol = 0;
		else f->vol += signum
				*sqrt(sqq)
				*t->fgs->vol
				/(cdim-depth+1);
	}
	vols(f, t->left, n, depth);
	vols(f, t->right, n, depth);

	return;
}


void find_volumes(fg *faces_gr, FILE *F) {
	if (!faces_gr) return;
	vols(faces_gr, faces_gr->facets, 0, 1);
	print_fg(faces_gr, F);
}


gsitef *get_site;
site_n *site_num;


void set_ch_root(simplex *s) {ch_root = s; return;}
/* set root to s, for purposes of getting normals etc. */


simplex *build_convex_hull(gsitef *get_s, site_n *site_numm, short dim, short vdd) {

/*
	get_s		returns next site each call;
			hull construction stops when NULL returned;
	site_numm	returns number of site when given site;
	dim		dimension of point set;
	vdd		if (vdd) then return Delaunay triangulation


*/

	simplex *s, *root;

	if (!Huge) Huge = DBL_MAX*DBL_MAX;

	cdim = 0;
	get_site = get_s;
	site_num = site_numm;
	pdim = dim;
	vd = vdd;

	exact_bits = DBL_MANT_DIG*log(FLT_RADIX)/log(2);
	b_err_min = DBL_EPSILON*MAXDIM*(1<<MAXDIM)*MAXDIM*3.01;
	b_err_min_sq = b_err_min * b_err_min;

	assert(get_site!=NULL); assert(site_num!=NULL);

	rdim = vd ? pdim+1 : pdim;
	if (rdim > MAXDIM)
		panic("dimension bound MAXDIM exceeded; rdim=%d; pdim=%d\n",
			rdim, pdim);
/*	fprintf(DFILE, "rdim=%d; pdim=%d\n", rdim, pdim); fflush(DFILE);*/

	point_size = site_size = sizeof(Coord)*pdim;
	basis_vec_size = sizeof(Coord)*rdim;
	basis_s_size = sizeof(basis_s)+ (2*rdim-1)*sizeof(Coord);
	simplex_size = sizeof(simplex) + (rdim-1)*sizeof(neighbor);
	Tree_size = sizeof(Tree);
	fg_size = sizeof(fg);

	root = NULL;
	if (vd) {
		p = infinity;
		NEWLRC(basis_s, infinity_basis);
		infinity_basis->vecs[2*rdim-1]
			= infinity_basis->vecs[rdim-1]
			= 1;
		infinity_basis->sqa
			= infinity_basis->sqb
			= 1;
	} else if (!(p = (*get_site)())) return 0;

	NEWL(simplex,root);

	ch_root = root;

	copy_simp(s,root);
	root->peak.vert = p;
	root->peak.simp = s;
	s->peak.simp = root;

	buildhull(root);
	return root;
}


void free_hull_storage(void) {
	free_basis_s_storage();
	free_simplex_storage();
	free_Tree_storage();
	free_fg_storage();
}


/* fg.c : face graph of hull, and splay trees */

/*
 * Ken Clarkson wrote this.  Copyright (c) 1995 by AT&T..
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */


#include <float.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include <string.h>

Tree* insert(site, Tree*);

Tree* tree_delete(site, Tree*);

void printtree(Tree*,int);

void printtree_flat(Tree*);


/* splay tree code */

/*
Fri Oct 21 21:15:01 EDT 1994
	style changes, removed Sedgewickized...
	Ken Clarkson

*/
/*
           An implementation of top-down splaying with sizes
             D. Sleator <sleator@cs.cmu.edu>, January 1994.

  This extends top-down-splay.c to maintain a size field in each node.
  This is the number of nodes in the subtree rooted there.  This makes
  it possible to efficiently compute the rank of a key.  (The rank is
  the number of nodes to the left of the given key.)  It it also
  possible to quickly find the node of a given rank.  Both of these
  operations are illustrated in the code below.  The remainder of this
  introduction is taken from top-down-splay.c.

  "Splay trees", or "self-adjusting search trees" are a simple and
  efficient data structure for storing an ordered set.  The data
  structure consists of a binary tree, with no additional fields.  It
  allows searching, insertion, deletion, deletemin, deletemax,
  splitting, joining, and many other operations, all with amortized
  logarithmic performance.  Since the trees adapt to the sequence of
  requests, their performance on real access patterns is typically even
  better.  Splay trees are described in a number of texts and papers
  [1,2,3,4].

  The code here is adapted from simple top-down splay, at the bottom of
  page 669 of [2].  It can be obtained via anonymous ftp from
  spade.pc.cs.cmu.edu in directory /usr/sleator/public.

  The chief modification here is that the splay operation works even if the
  item being splayed is not in the tree, and even if the tree root of the
  tree is NULL.  So the line:

                              t = splay(i, t);

  causes it to search for item with key i in the tree rooted at t.  If it's
  there, it is splayed to the root.  If it isn't there, then the node put
  at the root is the last one before NULL that would have been reached in a
  normal binary search for i.  (It's a neighbor of i in the tree.)  This
  allows many other operations to be easily implemented, as shown below.

  [1] "Data Structures and Their Algorithms", Lewis and Denenberg,
       Harper Collins, 1991, pp 243-251.
  [2] "Self-adjusting Binary Search Trees" Sleator and Tarjan,
       JACM Volume 32, No 3, July 1985, pp 652-686.
  [3] "Data Structure and Algorithm Analysis", Mark Weiss,
       Benjamin Cummins, 1992, pp 119-130.
  [4] "Data Structures, Algorithms, and Performance", Derick Wood,
       Addison-Wesley, 1993, pp 367-375
*/


STORAGE(Tree)


#define compare(i,j) (site_num(i)-site_num(j))
/* This is the comparison.                                       */
/* Returns <0 if i<j, =0 if i=j, and >0 if i>j                   */
 
#define node_size(x) ((x) ? ((x)->size) : 0 )
/* This macro returns the size of a node.  Unlike "x->size",     */
/* it works even if x=NULL.  The test could be avoided by using  */
/* a special version of NULL which was a real node with size 0.  */
 
Tree * splay (site i, Tree *t) 
/* Splay using the key i (which may or may not be in the tree.) */
/* The starting root is t, and the tree used is defined by rat  */
/* size fields are maintained */
{
    Tree N, *l, *r, *y;
    int comp, root_size, l_size, r_size;
    
    if (!t) return t;
    N.left = N.right = NULL;
    l = r = &N;
    root_size = node_size(t);
    l_size = r_size = 0;
 
    for (;;) {
        comp = compare(i, t->key);
        if (comp < 0) {
            if (!t->left) break;
            if (compare(i, t->left->key) < 0) {
                y = t->left;                           /* rotate right */
                t->left = y->right;
                y->right = t;
                t->size = node_size(t->left) + node_size(t->right) + 1;
                t = y;
                if (!t->left) break;
            }
            r->left = t;                               /* link right */
            r = t;
            t = t->left;
            r_size += 1+node_size(r->right);
        } else if (comp > 0) {
            if (!t->right) break;
            if (compare(i, t->right->key) > 0) {
                y = t->right;                          /* rotate left */
                t->right = y->left;
                y->left = t;
		t->size = node_size(t->left) + node_size(t->right) + 1;
                t = y;
                if (!t->right) break;
            }
            l->right = t;                              /* link left */
            l = t;
            t = t->right;
            l_size += 1+node_size(l->left);
        } else break;
    }
    l_size += node_size(t->left);  /* Now l_size and r_size are the sizes of */
    r_size += node_size(t->right); /* the left and right trees we just built.*/
    t->size = l_size + r_size + 1;

    l->right = r->left = NULL;

    /* The following two loops correct the size fields of the right path  */
    /* from the left child of the root and the right path from the left   */
    /* child of the root.                                                 */
    for (y = N.right; y != NULL; y = y->right) {
        y->size = l_size;
        l_size -= 1+node_size(y->left);
    }
    for (y = N.left; y != NULL; y = y->left) {
        y->size = r_size;
        r_size -= 1+node_size(y->right);
    }
 
    l->right = t->left;                                /* assemble */
    r->left = t->right;
    t->left = N.right;
    t->right = N.left;

    return t;
}

Tree * insert(site i, Tree * t) {
/* Insert key i into the tree t, if it is not already there. */
/* Return a pointer to the resulting tree.                   */
    Tree * new_site;

    if (t != NULL) {
	t = splay(i,t);
	if (compare(i, t->key)==0) {
	    return t;  /* it's already there */
	}
    }
    NEWL(Tree,new_site)
    if (!t) {
	new_site->left = new_site->right = NULL;
    } else if (compare(i, t->key) < 0) {
	new_site->left = t->left;
	new_site->right = t;
	t->left = NULL;
	t->size = 1+node_size(t->right);
    } else {
	new_site->right = t->right;
	new_site->left = t;
	t->right = NULL;
	t->size = 1+node_size(t->left);
    }
    new_site->key = i;
    new_site->size = 1 + node_size(new_site->left) + node_size(new_site->right);
    return new_site;
}

Tree * tree_delete(site i, Tree *t) {
/* Deletes i from the tree if it's there.               */
/* Return a pointer to the resulting tree.              */
    Tree * x;
    int tsize;

    if (!t) return NULL;
    tsize = t->size;
    t = splay(i,t);
    if (compare(i, t->key) == 0) {               /* found it */
	if (!t->left) {
	    x = t->right;
	} else {
	    x = splay(i, t->left);
	    x->right = t->right;
	}
	FREEL(Tree,t);
	if (x) x->size = tsize-1;
	return x;
    } else {
	return t;                         /* It wasn't there */
    }
}

Tree *find_rank(int r, Tree *t) {
/* Returns a pointer to the node in the tree with the given rank.  */
/* Returns NULL if there is no such node.                          */
/* Does not change the tree.  To guarantee logarithmic behavior,   */
/* the node found here should be splayed to the root.              */
    int lsize;
    if ((r < 0) || (r >= node_size(t))) return NULL;
    for (;;) {
	lsize = node_size(t->left);
	if (r < lsize) {
	    t = t->left;
	} else if (r > lsize) {
	    r = r - lsize -1;
	    t = t->right;
	} else {
	    return t;
	}
    }
}

void printtree_flat_inner(Tree * t) {
    if (!t) return;

    printtree_flat_inner(t->right);
    printf("%d ", t->key);fflush(stdout);
    printtree_flat_inner(t->left);
}

void printtree_flat(Tree * t) {
	if (!t) {
		printf("<empty tree>");
		return;
	}
	printtree_flat_inner(t);
}


void printtree(Tree * t, int d) {
    int i;
    if (!t) return;

   printtree(t->right, d+1);
    for (i=0; i<d; i++) printf("  ");
    printf("%p(%d)\n", t->key, t->size);fflush(stdout);
    printtree(t->left, d+1);
}






fg *faces_gr_t;

STORAGE(fg)

#define snkey(x) site_num((x)->vert)

fg *find_fg(simplex *s,int q) {

	fg *f;
	neighbor *si, *sn = s->neigh;
	Tree *t;

	if (q==0) return faces_gr_t;
	if (!faces_gr_t) NEWLRC(fg, faces_gr_t);
	f = faces_gr_t;
	for (si=sn; si<sn+cdim; si++) if (q & (1<<(si-sn))) {
		t = f->facets = insert(si->vert,f->facets);
		if (!t->fgs) NEWLRC(fg, (t->fgs))
		f = t->fgs;
	}
	return f;
}

void *add_to_fg(simplex *s, void *dum) {

	neighbor t, *si, *sj, *sn = s->neigh;
	fg *fq;
	int q,m,Q=1<<cdim, s1;
					/* sort neigh by site number */
	for (si=sn+2;si<sn+cdim;si++)
		for (sj=si; sj>sn+1 && snkey(sj-1) > snkey(sj); sj--)
			{t=*(sj-1); *(sj-1) = *sj; *sj = t;}

	NULLIFY(basis_s,s->normal);
	NULLIFY(basis_s,s->neigh[0].basis);

					/* insert subsets */
	for (q=1; q<Q; q++) find_fg(s,q);

					/* include all superset relations */
	for (q=1; q<Q; q++) {
		fq = find_fg(s,q);
		assert(fq);
		for (m=1,si=sn;si<sn+cdim;si++,m<<=1) if (!(q&m)) {
			fq->facets = insert(si->vert,fq->facets);
			fq->facets->fgs = find_fg(s, q|m);
		}
	}
	return NULL;	
}

fg *build_fg(simplex *root) {
	faces_gr_t= 0;
	visit_hull(root, add_to_fg);
	return faces_gr_t;
}

void visit_fg_i(   void (*v_fg)(Tree *, int, int),
			Tree *t, int depth, int vn, int boundary) {
	int	boundaryc = boundary;

	if (!t) return;

	assert(t->fgs);
	if (t->fgs->mark!=vn) {
		t->fgs->mark = vn;
		if (t->key!=infinity && !mo[site_num(t->key)]) boundaryc = 0; 
		v_fg(t,depth, boundaryc);
		visit_fg_i(v_fg, t->fgs->facets,depth+1, vn, boundaryc);
	}
	visit_fg_i(v_fg, t->left,depth,vn, boundary);
	visit_fg_i(v_fg, t->right,depth,vn,boundary);
}

void visit_fg(fg *faces_gr, void (*v_fg)(Tree *, int, int)) {
	static int fg_vn;
	visit_fg_i(v_fg, faces_gr->facets, 0, ++fg_vn, 1);
}

int visit_fg_i_far(void (*v_fg)(Tree *, int),
			Tree *t, int depth, int vn) {
	int nb = 0;

	if (!t) return 0;

	assert(t->fgs);
	if (t->fgs->mark!=vn) {
		t->fgs->mark = vn;
		nb = (t->key==infinity) || mo[site_num(t->key)];
		if (!nb && !visit_fg_i_far(v_fg, t->fgs->facets,depth+1,vn))
			v_fg(t,depth);
	}
	nb = visit_fg_i_far(v_fg, t->left,depth,vn) || nb;
	nb = visit_fg_i_far(v_fg, t->right,depth,vn) || nb;
	return nb;
}

void visit_fg_far(fg *faces_gr, void (*v_fg)(Tree *, int)) {
	static int fg_vn;
	visit_fg_i_far(v_fg,faces_gr->facets, 0, --fg_vn);
}



FILE *FG_OUT;

void p_fg(Tree* t, int depth, int bad) {
	static int fa[MAXDIM];
	int i;
	static double mults[MAXDIM];

	if (mults[0]==0) {
		mults[pdim] = 1;
		for (i=pdim-1; i>=0; i--) mults[i] = mult_up*mults[i+1];
	}

	fa[depth] = site_num(t->key);
	for (i=0;i<=depth;i++)
		fprintf(FG_OUT, "%d ", fa[i]);
	fprintf(FG_OUT, "	%G\n", t->fgs->vol/mults[depth]);
}

int p_fg_x_depth;

void p_fg_x(Tree*t, int depth, int bad) {

	static int fa[MAXDIM];
	static point fp[MAXDIM];
	int i,th;
	point tp;

	fa[depth] = site_num(t->key);
	fp[depth] = t->key;

	if (depth==p_fg_x_depth) for (i=0;i<=depth;i++)
			fprintf(FG_OUT, "%d%s", fa[i], (i==depth) ? "\n" : " ");
}

void print_fg_alt(fg *faces_gr, FILE *F, int fd) {
	FG_OUT=F;
 	if (!faces_gr) return;
	p_fg_x_depth = fd;
	visit_fg(faces_gr, p_fg_x);
	fclose(FG_OUT);
}


void print_fg(fg *faces_gr, FILE *F) {FG_OUT=F; visit_fg(faces_gr, p_fg);}


double fg_hist[100][100], fg_hist_bad[100][100],fg_hist_far[100][100];

void h_fg(Tree *t, int depth, int bad) {
	int i;
	if (!t->fgs->facets) return;
	if (bad) {
		fg_hist_bad[depth][t->fgs->facets->size]++;
		return;
	}
	fg_hist[depth][t->fgs->facets->size]++;
}

void h_fg_far(Tree* t, int depth) {
	if (t->fgs->facets) fg_hist_far[depth][t->fgs->facets->size]++;
}


void print_hist_fg(simplex *root, fg *faces_gr, FILE *F) {
	int i,j,k;
	double tot_good[100], tot_bad[100], tot_far[100];
	for (i=0;i<20;i++) {
		tot_good[i] = tot_bad[i] = tot_far[i] = 0;
		for (j=0;j<100;j++) {
			fg_hist[i][j]= fg_hist_bad[i][j]= fg_hist_far[i][j] = 0;
		}
	}
	if (!root) return;

	find_alpha(root);

	if (!faces_gr) faces_gr = build_fg(root);

	visit_fg(faces_gr, h_fg);
	visit_fg_far(faces_gr, h_fg_far);

	for (j=0;j<100;j++) for (i=0;i<20;i++) {
		tot_good[i] += fg_hist[i][j];
		tot_bad[i] += fg_hist_bad[i][j];
		tot_far[i]  += fg_hist_far[i][j];
	}

	for (i=19;i>=0 && !tot_good[i] && !tot_bad[i]; i--);
	fprintf(F,"totals	");
	for (k=0;k<=i;k++) {
		if (k==0) fprintf(F, "	");
		else fprintf(F,"			");
		fprintf(F, "%d/%d/%d",
			(int)tot_far[k], (int)tot_good[k], (int)tot_good[k] + (int)tot_bad[k]);
	}
		
	
	for (j=0;j<100;j++) {
		for (i=19; i>=0 && !fg_hist[i][j] && !fg_hist_bad[i][j]; i--);
		if (i==-1) continue;
		fprintf(F, "\n%d	",j);fflush(F);
			
		for (k=0;k<=i;k++) {
			if (k==0) fprintf(F, "	");
			else fprintf(F,"			");
			if (fg_hist[k][j] || fg_hist_bad[k][j])
				fprintf(F,
					"%2.1f/%2.1f/%2.1f",
					tot_far[k] ? 100*fg_hist_far[k][j]/tot_far[k]+.05 : 0,
					tot_good[k] ? 100*fg_hist[k][j]/tot_good[k]+.05 : 0,
					100*(fg_hist[k][j]+fg_hist_bad[k][j])/(tot_good[k]+tot_bad[k])+.05
				);
		}
	}
	fprintf(F, "\n");
}

 
/* io.c : input-output */


/*
 * Ken Clarkson wrote this.  Copyright (c) 1995 by AT&T..
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <float.h>
#include <math.h>

double mult_up = 1.0;
Coord mins[MAXDIM]
	= {DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX},
	maxs[MAXDIM]
	= {-DBL_MAX,-DBL_MAX,-DBL_MAX,-DBL_MAX,-DBL_MAX,-DBL_MAX,-DBL_MAX,-DBL_MAX};

void panic(char *fmt, ...) {
	va_list args;

	va_start(args, fmt);
	vfprintf(DFILE, fmt, args);
	fflush(DFILE);
	va_end(args);

	exit(1);
}

/*
FILE * popen(char*, char*);
void pclose(FILE*);
*/

char tmpfilenam[L_tmpnam];

FILE* efopen(char *file, char *mode) {
	FILE* fp;
	if (fp = fopen(file, mode)) return fp;
	fprintf(DFILE, "couldn't open file %s mode %s\n",file,mode);
	exit(1);
	return NULL;
}

FILE* epopen(char *com, char *mode) {
	FILE* fp;
	if (fp = popen(com, mode)) return fp;
	fprintf(stderr, "couldn't open stream %s mode %s\n",com,mode);
	exit(1);
	return 0;
}



void print_neighbor_snum(FILE* F, neighbor *n){

	assert(site_num!=NULL);
	if (n->vert)
		fprintf(F, "%d ", (*site_num)(n->vert));
	else
		fprintf(F, "NULL vert ");
	fflush(stdout);
}

void print_basis(FILE *F, basis_s* b) {
	if (!b) {fprintf(F, "NULL basis ");fflush(stdout);return;}
	if (b->lscale<0) {fprintf(F, "\nbasis computed");return;}
	fprintf(F, "\n%p  %d \n b=",(void*)b,b->lscale);
	print_point(F, rdim, b->vecs);
	fprintf(F, "\n a= ");
	print_point_int(F, rdim, b->vecs+rdim); fprintf(F, "   ");
	fflush(F);
}

void print_simplex_num(FILE *F, simplex *s) {
	fprintf(F, "simplex ");
	if(!s) fprintf(F, "NULL ");
	else fprintf(F, "%p  ", (void*)s);
}

void print_neighbor_full(FILE *F, neighbor *n){

	if (!n) {fprintf(F, "null neighbor\n");return;}

	print_simplex_num(F, n->simp);
	print_neighbor_snum(F, n);fprintf(F, ":  ");fflush(F);
	if (n->vert) {
/*		if (n->basis && n->basis->lscale <0) fprintf(F, "trans ");*/
		/* else */ print_point(F, pdim,n->vert);
		fflush(F);
	}
	print_basis(F, n->basis);
	fflush(F);
	fprintf(F, "\n");
}

void *print_facet(FILE *F, simplex *s, print_neighbor_f *pnfin) {
	int i;
	neighbor *sn = s->neigh;

/*	fprintf(F, "%d ", s->mark);*/
	for (i=0; i<cdim;i++,sn++) (*pnfin)(F, sn);
	fprintf(F, "\n");
	fflush(F);
	return NULL;
}




void *print_simplex_f(simplex *s, FILE *F, print_neighbor_f *pnfin){

	static print_neighbor_f *pnf;

	if (pnfin) {pnf=pnfin; if (!s) return NULL;}

	print_simplex_num(F, s);
	fprintf(F, "\n");
	if (!s) return NULL;
	fprintf(F, "normal ="); print_basis(F, s->normal); fprintf(F, "\n");
	fprintf(F, "peak ="); (*pnf)(F, &(s->peak));
	fprintf (F, "facet =\n");fflush(F);
	return print_facet(F, s, pnf);
}

void *print_simplex(simplex *s, void *Fin) {
	static FILE *F;

	if (Fin) {F=(FILE*)Fin; if (!s) return NULL;}

	return print_simplex_f(s, F, 0);

}


void print_triang(simplex *root, FILE *F, print_neighbor_f *pnf) {
	print_simplex(0,F);
	print_simplex_f(0,0,pnf);
	visit_triang(root, print_simplex);
}

void *p_peak_test(simplex *s) {return (s->peak.vert==p) ? (void*)s : (void*)NULL;}



void *check_simplex(simplex *s, void *dum){

	int i,j,k,l;
	neighbor *sn, *snn, *sn2;
	simplex *sns;
	site vn;

	for (i=-1,sn=s->neigh-1;i<cdim;i++,sn++) {
		sns = sn->simp;
		if (!sns) {
			fprintf(DFILE, "check_triang; bad simplex\n");
			print_simplex_f(s, DFILE, &print_neighbor_full); fprintf(DFILE, "site_num(p)=%G\n",site_num(p));
			return s;
		}
		if (!s->peak.vert && sns->peak.vert && i!=-1) {
			fprintf(DFILE, "huh?\n");
			print_simplex_f(s, DFILE, &print_neighbor_full);
			print_simplex_f(sns, DFILE, &print_neighbor_full);
			exit(1);
		}
		for (j=-1,snn=sns->neigh-1; j<cdim && snn->simp!=s; j++,snn++);
		if (j==cdim) {
			fprintf(DFILE, "adjacency failure:\n");
			DEBEXP(-1,site_num(p))
			print_simplex_f(sns, DFILE, &print_neighbor_full);
			print_simplex_f(s, DFILE, &print_neighbor_full);
			exit(1);
		}
		for (k=-1,snn=sns->neigh-1; k<cdim; k++,snn++){
			vn = snn->vert;
			if (k!=j) {
				for (l=-1,sn2 = s->neigh-1;
					l<cdim && sn2->vert != vn;
					l++,sn2++);
				if (l==cdim) {
					fprintf(DFILE, "cdim=%d\n",cdim);
					fprintf(DFILE, "error: neighboring simplices with incompatible vertices:\n");
					print_simplex_f(sns, DFILE, &print_neighbor_full);
					print_simplex_f(s, DFILE, &print_neighbor_full);
					exit(1);
				}	
			}
		}
	}
	return NULL;
}

int p_neight(simplex *s, int i, void *dum) {return s->neigh[i].vert !=p;}

void check_triang(simplex *root){visit_triang(root, &check_simplex);}

void check_new_triangs(simplex *s){visit_triang_gen(s, check_simplex, p_neight);}





/* outfuncs: given a list of points, output in a given format */


void vlist_out(point *v, int vdim, FILE *Fin, int amble) {

	static FILE *F;
	int j;

	if (Fin) {F=Fin; if (!v) return;}

	for (j=0;j<vdim;j++) fprintf(F, "%d ", (site_num)(v[j]));
	fprintf(F,"\n");

	return;
}

void off_out(point *v, int vdim, FILE *Fin, int amble) {

	static FILE *F, *G;
	static FILE *OFFFILE;
	static char offfilenam[L_tmpnam];
	char comst[100], buf[200];
	int j,i;

	if (Fin) {F=Fin;}

	if (pdim!=3) { warning(-10, off apparently for 3d points only); return;}

	if (amble==0) {
		for (i=0;i<vdim;i++) if (v[i]==infinity) return;
		fprintf(OFFFILE, "%d ", vdim);
		for (j=0;j<vdim;j++) fprintf(OFFFILE, "%d ", (site_num)(v[j]));
		fprintf(OFFFILE,"\n");
	} else if (amble==-1) {
		OFFFILE = efopen(tmpnam(offfilenam), "w");
	} else {
		fclose(OFFFILE);

		fprintf(F, "	OFF\n");
	
		sprintf(comst, "wc %s", tmpfilenam);
		G = epopen(comst, "r");
		fscanf(G, "%d", &i);
		fprintf(F, " %d", i);
		pclose(G);
	
		sprintf(comst, "wc %s", offfilenam);
		G = epopen(comst, "r");
		fscanf(G, "%d", &i);
		fprintf(F, " %d", i);
		pclose(G);
	
		fprintf (F, " 0\n");
	
		G = efopen(tmpfilenam, "r");
		while (fgets(buf, sizeof(buf), G)) fprintf(F, "%s", buf);
		fclose(G);
	
		G = efopen(offfilenam, "r");
	
	
		while (fgets(buf, sizeof(buf), G)) fprintf(F, "%s", buf);
		fclose(G);
	}

}



void mp_out(point *v, int vdim, FILE *Fin, int amble) {


/* should fix scaling */

	static int figno=1;
	static FILE *F;

	if (Fin) {F=Fin;}

	if (pdim!=2) { warning(-10, mp for planar points only); return;}
	if (amble==0) {
		int i;
		if (!v) return;
		for (i=0;i<vdim;i++) if (v[i]==infinity) {
			point t=v[i];
			v[i]=v[vdim-1];
			v[vdim-1] = t;
			vdim--;
			break;
		}
		fprintf(F, "draw ");
		for (i=0;i<vdim;i++) 
			fprintf(F,
				(i+1<vdim) ? "(%Gu,%Gu)--" : "(%Gu,%Gu);\n",
				v[i][0]/mult_up,v[i][1]/mult_up
			);
	} else if (amble==-1) {
		if (figno==1) fprintf(F, "u=1pt;\n");
		fprintf(F , "beginfig(%d);\n",figno++);
	} else if (amble==1) {
		fprintf(F , "endfig;\n");
	}
}


void ps_out(point *v, int vdim, FILE *Fin, int amble) {

	static FILE *F;
	static double scaler;

	if (Fin) {F=Fin;}

	if (pdim!=2) { warning(-10, ps for planar points only); return;}

	if (amble==0) {
		int i;
		if (!v) return;
		for (i=0;i<vdim;i++) if (v[i]==infinity) {
			point t=v[i];
			v[i]=v[vdim-1];
			v[vdim-1] = t;
			vdim--;
			break;
		}
		fprintf(F,
			"newpath %G %G moveto\n",
			v[0][0]*scaler,v[0][1]*scaler);
		for (i=1;i<vdim;i++) 
			fprintf(F,
				"%G %G lineto\n",
				v[i][0]*scaler,v[i][1]*scaler
			);
		fprintf(F, "stroke\n");
	} else if (amble==-1) {
		float len[2], maxlen;
		fprintf(F, "%%!PS\n");
		len[0] = maxs[0]-mins[0]; len[1] = maxs[1]-mins[1];
		maxlen = (len[0]>len[1]) ? len[0] : len[1];
		scaler = 216/maxlen;
	
		fprintf(F, "%%%%BoundingBox: %G %G %G %G \n",
			mins[0]*scaler,
			mins[1]*scaler,
			maxs[0]*scaler,
			maxs[1]*scaler);
		fprintf(F, "%%%%Creator: hull program\n");
		fprintf(F, "%%%%Pages: 1\n");
		fprintf(F, "%%%%EndProlog\n");
		fprintf(F, "%%%%Page: 1 1\n");
		fprintf(F, " 0.5 setlinewidth [] 0 setdash\n");
		fprintf(F, " 1 setlinecap 1 setlinejoin 10 setmiterlimit\n");
	} else if (amble==1) {
		fprintf(F , "showpage\n %%%%EOF\n");
	}
}

void cpr_out(point *v, int vdim, FILE *Fin, int amble) {

	static FILE *F;
	int i;

	if (Fin) {F=Fin; if (!v) return;}

	if (pdim!=3) { warning(-10, cpr for 3d points only); return;}
	
	for (i=0;i<vdim;i++) if (v[i]==infinity) return;

	fprintf(F, "t %G %G %G %G %G %G %G %G %G 3 128\n",
		v[0][0]/mult_up,v[0][1]/mult_up,v[0][2]/mult_up,
		v[1][0]/mult_up,v[1][1]/mult_up,v[1][2]/mult_up,
		v[2][0]/mult_up,v[2][1]/mult_up,v[2][2]/mult_up
	);
}


/* vist_funcs for different kinds of output: facets, alpha shapes, etc. */





void *facets_print(simplex *s, void *p) {

	static out_func *out_func_here;
	point v[MAXDIM];
	int j;

	if (p) {out_func_here = (out_func*)p; if (!s) return NULL;}

	for (j=0;j<cdim;j++) v[j] = s->neigh[j].vert;

	out_func_here(v,cdim,0,0);

	return NULL;
}


void *ridges_print(simplex *s, void *p) {

	static out_func *out_func_here;
	point v[MAXDIM];
	int j,k,vnum;

	if (p) {out_func_here = (out_func*)p; if (!s) return NULL;}

	for (j=0;j<cdim;j++) {
		vnum=0;
		for (k=0;k<cdim;k++) {
			if (k==j) continue;
			v[vnum++] = (s->neigh[k].vert);
		}
		out_func_here(v,cdim-1,0,0);
	}
	return NULL;
}



void *afacets_print(simplex *s, void *p) {

	static out_func *out_func_here;
	point v[MAXDIM];
	int j,k,vnum;

	if (p) {out_func_here = (out_func*)p; if (!s) return NULL;}

	for (j=0;j<cdim;j++) { /* check for ashape consistency */
		for (k=0;k<cdim;k++) if (s->neigh[j].simp->neigh[k].simp==s) break;
		if (alph_test(s,j,0)!=alph_test(s->neigh[j].simp,k,0)) {
			DEB(-10,alpha-shape not consistent)
			DEBTR(-10)
			print_simplex_f(s,DFILE,&print_neighbor_full);
			print_simplex_f(s->neigh[j].simp,DFILE,&print_neighbor_full);
			fflush(DFILE);
			exit(1);
		}
	}
	for (j=0;j<cdim;j++) {
		vnum=0;
		if (alph_test(s,j,0)) continue;
		for (k=0;k<cdim;k++) {
			if (k==j) continue;
			v[vnum++] = s->neigh[k].vert;
		}
		out_func_here(v,cdim-1,0,0);
	}
	return NULL;
}


/* pointops.c */

/*
 * Ken Clarkson wrote this.  Copyright (c) 1995 by AT&T..
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */

#include <float.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

int	pdim;	/* point dimension */

#define NEARZERO(d)	((d) < FLT_EPSILON && (d) > -FLT_EPSILON)
Coord maxdist(int dim, point p1, point p2) {
	Coord	x,y,
		d = 0;
	int i = dim;


	while (i--) {
		x = *p1++;
		y = *p2++;
		d += (x<y) ? y-x : x-y ;
	}

	return d;
}

void print_point(FILE *F, int dim, point p) {
	int j;
	if (!p) {
		fprintf(F, "NULL");
		return;
	}
	for (j=0;j<dim;j++) fprintf(F, "%g  ", *p++);
}

void print_point_int(FILE *F, int dim, point p) {
	int j;
	if (!p) {
		fprintf(F, "NULL");
		return;
	}
	for (j=0;j<dim;j++) fprintf(F, "%.20g  ", *p++);
}


int scale(int dim, point p) {
	Coord max = 0;
	int i;
	Coord abs,val;
	for (i=0;i<dim;i++) {
		val = p[i];
		abs = (val > 0) ? val: -val;
		max = (abs > max) ? abs : max;
	}


	if (max< 100*DBL_EPSILON) {
		fprintf(stderr, "fails to scale: ");
		print_point(stderr, dim,p);fflush(stderr);
		fprintf(stderr, "\n");
		return 1;
	}

	for (i=0;i<dim;i++) p[i] /= max;

	return 0;
}


/*
int normalize(int dim, point p) {
	Coord norm;
	int i;

	norm = norm2(dim,p);

	if (norm < 3*FLT_EPSILON) {
		fprintf(stderr, "fails to normalize: ");
		print_point(dim,p);fflush(stdout);
		fprintf(stderr, "\n");
		return 1;
	}


	for (i=0;i<dim;i++) p[i] /= norm;

	return 0;
}

*/

/* clarkson_convex_hull */
/* by Rephael Wenger, November, 2001 */

FILE * DFILE;

static int cch_dimension = 0;
static Coord * cch_point_coord = NULL;
static int cch_num_points = 0;
static int cch_current_point = 0;

static int cch_convex_hull_dimension = 0;
static int cch_num_simplices = 0;
static int * cch_current_simplex_vert = NULL;
static int * cch_simplex_vert = NULL;

site clarkson_convex_hull_next_site(void);
long clarkson_convex_hull_site_number(site p);
void * clarkson_convex_hull_count_simplex(simplex *s, void * ignore);
void * clarkson_convex_hull_add_simplex(simplex *s, void * ignore);


void clarkson_convex_hull
  (int dimension, double * point_coord, int num_points, 
   int * convex_hull_dimension, int ** simplex_vert, int * num_simplices)
// point_coord = list of point coordinates
// num_points = number of points
// convex_hull_dimension = dimension of convex hull (including interior)
// simplex_vert = list of simplex vertices
// num_simplices = number of simplices
{
  struct simplex * root;
  int num_simplex_vert;            // total number of simplex vertices

  DFILE = stderr;

  cch_dimension = dimension;
  cch_point_coord = point_coord;
  cch_num_points = num_points;

  (*convex_hull_dimension) = 0;
  (*simplex_vert) = NULL;
  (*num_simplices) = 0;

  if (dimension < 1 || num_points == 0)
    return;

  assert(point_coord != NULL);

  cch_current_point = 0;
  root = build_convex_hull
    (clarkson_convex_hull_next_site, clarkson_convex_hull_site_number,
     cch_dimension, 0);

  // initialize simplex data
  cch_convex_hull_dimension = cdim;
  cch_num_simplices = 0;
  cch_simplex_vert = NULL;

  // count number of simplices
  visit_hull(root, clarkson_convex_hull_count_simplex);

  // allocate memory
  num_simplex_vert = cch_num_simplices * cch_convex_hull_dimension;

  cch_simplex_vert = (int *) malloc(sizeof(int) * num_simplex_vert);

  // store simplex vertices
  cch_current_simplex_vert = cch_simplex_vert;
  visit_hull(root, clarkson_convex_hull_add_simplex);

  (*convex_hull_dimension) = cch_convex_hull_dimension;
  (*simplex_vert) = cch_simplex_vert;
  (*num_simplices) = cch_num_simplices;

  free_hull_storage();
}

site clarkson_convex_hull_next_site(void) 
// return next site
{
  Coord * point_coord_ptr;
  
  if (cch_current_point >= cch_num_points)
    return(NULL);

  point_coord_ptr =
    cch_point_coord + cch_current_point*cch_dimension;
  cch_current_point++;

  return(point_coord_ptr);
}

long clarkson_convex_hull_site_number(site p)
// return site number of p
{
  int n = (p-cch_point_coord)/cch_dimension;
  return(n);
}

void * clarkson_convex_hull_add_simplex(simplex *s, void * ignore)
// s = simplex
// ignore = ignore this field
{
  int j;
  for (j=0; j < cch_convex_hull_dimension; j++) {
    assert(cch_current_simplex_vert-cch_simplex_vert < 
	   cch_num_simplices*cch_convex_hull_dimension);

    *(cch_current_simplex_vert) = 
      clarkson_convex_hull_site_number(s->neigh[j].vert);
    cch_current_simplex_vert++;
  };

  return NULL;
}

void * clarkson_convex_hull_count_simplex(simplex *s, void * ignore)
// s = simplex
// ignore = ignore this field
{
  cch_num_simplices++;

  return NULL;
}

void free_simplex_vertices(int * simplex_vert)
{
  free(simplex_vert);
}

