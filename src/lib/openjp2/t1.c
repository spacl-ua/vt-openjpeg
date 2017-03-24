/*
 * The copyright in this software is being made available under the 2-clauses
 * BSD License, included below. This software may be subject to other third
 * party and contributor rights, including patent rights, and no such rights
 * are granted under this license.
 * Copyright (c) 2015-2017, The University of Arizona, U.S.A.
 * Copyright (c) 2002-2014, Universite catholique de Louvain (UCL), Belgium
 * Copyright (c) 2002-2014, Professor Benoit Macq
 * Copyright (c) 2001-2003, David Janssens
 * Copyright (c) 2002-2003, Yannick Verschueren
 * Copyright (c) 2003-2007, Francois-Olivier Devaux
 * Copyright (c) 2003-2014, Antonin Descampe
 * Copyright (c) 2005, Herve Drolon, FreeImage Team
 * Copyright (c) 2007, Callum Lerwick <seg@haxxed.com>
 * 
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS `AS IS'
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include "opj_includes.h"
#include "t1_luts.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*#define IS_DEBUG_PRINTF*/
#define IS_USING_MASKING 0
#define TSI_ENC 0
/*
 * Normal cdf for
 * {sqrt(32) < x < GAUSS_XUPPER} union { GAUSS_XLOWER < x < -sqrt(32) }.
 */


/** @defgroup T1 T1 - Implementation of the tier-1 coding */
/*@{*/

/** @name Local static functions */
/*@{*/

//visually lossless coding

static void yl_opj_t1_enc_sigpass_step(opj_t1_t *t1,
                                       opj_flag_t *flagsp,
                                       OPJ_INT32 *datap,
                                       OPJ_INT32 *dec_datap,
                                       OPJ_UINT32 orient,
                                       OPJ_INT32 bpno,
                                       OPJ_INT32 one,
                                       OPJ_INT32 *nmsedec,
                                       OPJ_BYTE type,
                                       OPJ_UINT32 vsc);


static void yl_opj_t1_enc_sigpass( opj_t1_t *t1,
                                   OPJ_INT32 *dec_datap,
                                   OPJ_INT32 bpno,
                                   OPJ_UINT32 orient,
                                   OPJ_INT32 *nmsedec,
                                   OPJ_BYTE type,
                                   OPJ_UINT32 cblksty);

void yl_opj_t1_enc_refpass(
    opj_t1_t *t1,
    OPJ_INT32 *datap_dec,
    OPJ_INT32 bpno,
    OPJ_INT32 *nmsedec,
    OPJ_BYTE type,
    OPJ_UINT32 cblksty);

void yl_opj_t1_enc_refpass_step( opj_t1_t *t1,
                                 opj_flag_t *flagsp,
                                 OPJ_INT32 *datap,
                                 OPJ_INT32 *datap_dec,
                                 OPJ_INT32 bpno,
                                 OPJ_INT32 one,
                                 OPJ_INT32 *nmsedec,
                                 OPJ_BYTE type,
                                 OPJ_UINT32 vsc);

void yl_opj_t1_enc_clnpass(
    opj_t1_t *t1,
    OPJ_INT32 *datap_dec,
    OPJ_INT32 bpno,
    OPJ_UINT32 orient,
    OPJ_INT32 *nmsedec,
    OPJ_UINT32 cblksty);

static void yl_opj_t1_encode_cblk(opj_t1_t *t1,
                                  opj_tcd_cblk_enc_t* cblk,
                                  OPJ_UINT32 orient,
                                  OPJ_UINT32 compno,
                                  OPJ_UINT32 level,
                                  OPJ_UINT32 qmfbid,
                                  OPJ_FLOAT64 stepsize,
                                  OPJ_UINT32 cblksty,
                                  OPJ_UINT32 numcomps,
                                  opj_tcd_tile_t * tile,
                                  const OPJ_FLOAT64 * mct_norms,
                                  OPJ_FLOAT64 *lin_mod_slope,
                                  OPJ_FLOAT64 *lin_mod_intercept,
                                  OPJ_INT32 num_vts,
                                  OPJ_INT32* tiledp
                                 );
//visually lossless coding :end

static INLINE OPJ_BYTE opj_t1_getctxno_zc(OPJ_UINT32 f, OPJ_UINT32 orient);
static OPJ_BYTE opj_t1_getctxno_sc(OPJ_UINT32 f);
static INLINE OPJ_UINT32 opj_t1_getctxno_mag(OPJ_UINT32 f);
static OPJ_BYTE opj_t1_getspb(OPJ_UINT32 f);
static OPJ_INT16 opj_t1_getnmsedec_sig(OPJ_UINT32 x, OPJ_UINT32 bitpos);
static OPJ_INT16 opj_t1_getnmsedec_ref(OPJ_UINT32 x, OPJ_UINT32 bitpos);
static void opj_t1_updateflags(opj_flag_t *flagsp, OPJ_UINT32 s, OPJ_UINT32 stride);
/**
Encode significant pass
*/
static void opj_t1_enc_sigpass_step(opj_t1_t *t1,
                                    opj_flag_t *flagsp,
                                    OPJ_INT32 *datap,
                                    OPJ_UINT32 orient,
                                    OPJ_INT32 bpno,
                                    OPJ_INT32 one,
                                    OPJ_INT32 *nmsedec,
                                    OPJ_BYTE type,
                                    OPJ_UINT32 vsc);






/**
Decode significant pass
*/
#if 0
static void opj_t1_dec_sigpass_step(opj_t1_t *t1,
                                    opj_flag_t *flagsp,
                                    OPJ_INT32 *datap,
                                    OPJ_UINT32 orient,
                                    OPJ_INT32 oneplushalf,
                                    OPJ_BYTE type,
                                    OPJ_UINT32 vsc);
#endif

static INLINE void opj_t1_dec_sigpass_step_raw(
    opj_t1_t *t1,
    opj_flag_t *flagsp,
    OPJ_INT32 *datap,
    OPJ_INT32 orient,
    OPJ_INT32 oneplushalf,
    OPJ_INT32 vsc);
static INLINE void opj_t1_dec_sigpass_step_mqc(
    opj_t1_t *t1,
    opj_flag_t *flagsp,
    OPJ_INT32 *datap,
    OPJ_INT32 orient,
    OPJ_INT32 oneplushalf);
static INLINE void opj_t1_dec_sigpass_step_mqc_vsc(
    opj_t1_t *t1,
    opj_flag_t *flagsp,
    OPJ_INT32 *datap,
    OPJ_INT32 orient,
    OPJ_INT32 oneplushalf,
    OPJ_INT32 vsc);


/**
Encode significant pass
*/
static void opj_t1_enc_sigpass( opj_t1_t *t1,
                                OPJ_INT32 bpno,
                                OPJ_UINT32 orient,
                                OPJ_INT32 *nmsedec,
                                OPJ_BYTE type,
                                OPJ_UINT32 cblksty);

/**
Decode significant pass
*/
static void opj_t1_dec_sigpass_raw(
    opj_t1_t *t1,
    OPJ_INT32 bpno,
    OPJ_INT32 orient,
    OPJ_INT32 cblksty);
static void opj_t1_dec_sigpass_mqc(
    opj_t1_t *t1,
    OPJ_INT32 bpno,
    OPJ_INT32 orient);
static void opj_t1_dec_sigpass_mqc_vsc(
    opj_t1_t *t1,
    OPJ_INT32 bpno,
    OPJ_INT32 orient);



/**
Encode refinement pass
*/
static void opj_t1_enc_refpass_step(opj_t1_t *t1,
                                    opj_flag_t *flagsp,
                                    OPJ_INT32 *datap,
                                    OPJ_INT32 bpno,
                                    OPJ_INT32 one,
                                    OPJ_INT32 *nmsedec,
                                    OPJ_BYTE type,
                                    OPJ_UINT32 vsc);


/**
Encode refinement pass
*/
static void opj_t1_enc_refpass( opj_t1_t *t1,
                                OPJ_INT32 bpno,
                                OPJ_INT32 *nmsedec,
                                OPJ_BYTE type,
                                OPJ_UINT32 cblksty);

/**
Decode refinement pass
*/
static void opj_t1_dec_refpass_raw(
    opj_t1_t *t1,
    OPJ_INT32 bpno,
    OPJ_INT32 cblksty);
static void opj_t1_dec_refpass_mqc(
    opj_t1_t *t1,
    OPJ_INT32 bpno);
static void opj_t1_dec_refpass_mqc_vsc(
    opj_t1_t *t1,
    OPJ_INT32 bpno);


/**
Decode refinement pass
*/
#if 0
static void opj_t1_dec_refpass_step(opj_t1_t *t1,
                                    opj_flag_t *flagsp,
                                    OPJ_INT32 *datap,
                                    OPJ_INT32 poshalf,
                                    OPJ_INT32 neghalf,
                                    OPJ_BYTE type,
                                    OPJ_UINT32 vsc);
#endif

static INLINE void  opj_t1_dec_refpass_step_raw(
    opj_t1_t *t1,
    opj_flag_t *flagsp,
    OPJ_INT32 *datap,
    OPJ_INT32 poshalf,
    OPJ_INT32 neghalf,
    OPJ_INT32 vsc);
static INLINE void opj_t1_dec_refpass_step_mqc(
    opj_t1_t *t1,
    opj_flag_t *flagsp,
    OPJ_INT32 *datap,
    OPJ_INT32 poshalf,
    OPJ_INT32 neghalf);
static INLINE void opj_t1_dec_refpass_step_mqc_vsc(
    opj_t1_t *t1,
    opj_flag_t *flagsp,
    OPJ_INT32 *datap,
    OPJ_INT32 poshalf,
    OPJ_INT32 neghalf,
    OPJ_INT32 vsc);



/**
Encode clean-up pass
*/
static void opj_t1_enc_clnpass_step(
    opj_t1_t *t1,
    opj_flag_t *flagsp,
    OPJ_INT32 *datap,
    OPJ_UINT32 orient,
    OPJ_INT32 bpno,
    OPJ_INT32 one,
    OPJ_INT32 *nmsedec,
    OPJ_UINT32 partial,
    OPJ_UINT32 vsc);
/**
Decode clean-up pass
*/
static void opj_t1_dec_clnpass_step_partial(
    opj_t1_t *t1,
    opj_flag_t *flagsp,
    OPJ_INT32 *datap,
    OPJ_INT32 orient,
    OPJ_INT32 oneplushalf);
static void opj_t1_dec_clnpass_step(
    opj_t1_t *t1,
    opj_flag_t *flagsp,
    OPJ_INT32 *datap,
    OPJ_INT32 orient,
    OPJ_INT32 oneplushalf);
static void opj_t1_dec_clnpass_step_vsc(
    opj_t1_t *t1,
    opj_flag_t *flagsp,
    OPJ_INT32 *datap,
    OPJ_INT32 orient,
    OPJ_INT32 oneplushalf,
    OPJ_INT32 partial,
    OPJ_INT32 vsc);
/**
Encode clean-up pass
*/
static void opj_t1_enc_clnpass(
    opj_t1_t *t1,
    OPJ_INT32 bpno,
    OPJ_UINT32 orient,
    OPJ_INT32 *nmsedec,
    OPJ_UINT32 cblksty);
/**
Decode clean-up pass
*/
static void opj_t1_dec_clnpass(
    opj_t1_t *t1,
    OPJ_INT32 bpno,
    OPJ_INT32 orient,
    OPJ_INT32 cblksty);

static OPJ_FLOAT64 opj_t1_getwmsedec(
    OPJ_INT32 nmsedec,
    OPJ_UINT32 compno,
    OPJ_UINT32 level,
    OPJ_UINT32 orient,
    OPJ_INT32 bpno,
    OPJ_UINT32 qmfbid,
    OPJ_FLOAT64 stepsize,
    OPJ_UINT32 numcomps,
    const OPJ_FLOAT64 * mct_norms);



static void opj_t1_encode_cblk( opj_t1_t *t1,
                                opj_tcd_cblk_enc_t* cblk,
                                OPJ_UINT32 orient,
                                OPJ_UINT32 compno,
                                OPJ_UINT32 level,
                                OPJ_UINT32 qmfbid,
                                OPJ_FLOAT64 stepsize,
                                OPJ_UINT32 cblksty,
                                OPJ_UINT32 numcomps,
                                opj_tcd_tile_t * tile,
                                const OPJ_FLOAT64 * mct_norms);

/**
Decode 1 code-block
@param t1 T1 handle
@param cblk Code-block coding parameters
@param orient
@param roishift Region of interest shifting value
@param cblksty Code-block style
*/
static OPJ_BOOL opj_t1_decode_cblk( opj_t1_t *t1,
                                    opj_tcd_cblk_dec_t* cblk,
                                    OPJ_UINT32 orient,
                                    OPJ_UINT32 roishift,
                                    OPJ_UINT32 cblksty);

OPJ_BOOL opj_t1_allocate_buffers(   opj_t1_t *t1,
                                    OPJ_UINT32 w,
                                    OPJ_UINT32 h);

/*@}*/

/*@}*/

/* ----------------------------------------------------------------------- */

OPJ_BYTE opj_t1_getctxno_zc(OPJ_UINT32 f, OPJ_UINT32 orient)
{
    return lut_ctxno_zc[(orient << 8) | (f & T1_SIG_OTH)];
}

OPJ_BYTE opj_t1_getctxno_sc(OPJ_UINT32 f)
{
    return lut_ctxno_sc[(f & (T1_SIG_PRIM | T1_SGN)) >> 4];
}

OPJ_UINT32 opj_t1_getctxno_mag(OPJ_UINT32 f)
{
    OPJ_UINT32 tmp1 = (f & T1_SIG_OTH) ? T1_CTXNO_MAG + 1 : T1_CTXNO_MAG;
    OPJ_UINT32 tmp2 = (f & T1_REFINE) ? T1_CTXNO_MAG + 2 : tmp1;
    return (tmp2);
}

OPJ_BYTE opj_t1_getspb(OPJ_UINT32 f)
{
    return lut_spb[(f & (T1_SIG_PRIM | T1_SGN)) >> 4];
}

OPJ_INT16 opj_t1_getnmsedec_sig(OPJ_UINT32 x, OPJ_UINT32 bitpos)
{
    if (bitpos > T1_NMSEDEC_FRACBITS)
    {
        return lut_nmsedec_sig[(x >> (bitpos - T1_NMSEDEC_FRACBITS)) & ((1 << T1_NMSEDEC_BITS) - 1)];
    }

    return lut_nmsedec_sig0[x & ((1 << T1_NMSEDEC_BITS) - 1)];
}

OPJ_INT16 opj_t1_getnmsedec_ref(OPJ_UINT32 x, OPJ_UINT32 bitpos)
{
    if (bitpos > T1_NMSEDEC_FRACBITS)
    {
        return lut_nmsedec_ref[(x >> (bitpos - T1_NMSEDEC_FRACBITS)) & ((1 << T1_NMSEDEC_BITS) - 1)];
    }

    return lut_nmsedec_ref0[x & ((1 << T1_NMSEDEC_BITS) - 1)];
}

void opj_t1_updateflags(opj_flag_t *flagsp, OPJ_UINT32 s, OPJ_UINT32 stride)
{
    opj_flag_t *np = flagsp - stride;
    opj_flag_t *sp = flagsp + stride;

    static const opj_flag_t mod[] =
    {
        T1_SIG_S, T1_SIG_S|T1_SGN_S,
        T1_SIG_E, T1_SIG_E|T1_SGN_E,
        T1_SIG_W, T1_SIG_W|T1_SGN_W,
        T1_SIG_N, T1_SIG_N|T1_SGN_N
    };

    np[-1] |= T1_SIG_SE;
    np[0]  |= mod[s];
    np[1]  |= T1_SIG_SW;

    flagsp[-1] |= mod[s+2];
    flagsp[0]  |= T1_SIG;
    flagsp[1]  |= mod[s+4];

    sp[-1] |= T1_SIG_NE;
    sp[0]  |= mod[s+6];
    sp[1]  |= T1_SIG_NW;
}

void opj_t1_enc_sigpass_step(   opj_t1_t *t1,
                                opj_flag_t *flagsp,
                                OPJ_INT32 *datap,
                                OPJ_UINT32 orient,
                                OPJ_INT32 bpno,
                                OPJ_INT32 one,
                                OPJ_INT32 *nmsedec,
                                OPJ_BYTE type,
                                OPJ_UINT32 vsc
                            )
{
    OPJ_INT32 v;
    OPJ_UINT32 flag;

    opj_mqc_t *mqc = t1->mqc;	/* MQC component */

    flag = vsc ? (OPJ_UINT32)((*flagsp) & (~(T1_SIG_S | T1_SIG_SE | T1_SIG_SW | T1_SGN_S))) : (OPJ_UINT32)(*flagsp);
    if ((flag & T1_SIG_OTH) && !(flag & (T1_SIG | T1_VISIT)))
    {
        v = opj_int_abs(*datap) & one ? 1 : 0;
        opj_mqc_setcurctx(mqc, opj_t1_getctxno_zc(flag, orient));	/* ESSAI */
        if (type == T1_TYPE_RAW)  	/* BYPASS/LAZY MODE */
        {
            opj_mqc_bypass_enc(mqc, (OPJ_UINT32)v);
        }
        else
        {
            opj_mqc_encode(mqc, (OPJ_UINT32)v);
        }
        if (v)
        {
            v = *datap < 0 ? 1 : 0;
            *nmsedec +=	opj_t1_getnmsedec_sig((OPJ_UINT32)opj_int_abs(*datap), (OPJ_UINT32)(bpno + T1_NMSEDEC_FRACBITS));
            opj_mqc_setcurctx(mqc, opj_t1_getctxno_sc(flag));	/* ESSAI */
            if (type == T1_TYPE_RAW)  	/* BYPASS/LAZY MODE */
            {
                opj_mqc_bypass_enc(mqc, (OPJ_UINT32)v);
            }
            else
            {
                opj_mqc_encode(mqc, (OPJ_UINT32)(v ^ opj_t1_getspb((OPJ_UINT32)flag)));
            }
            opj_t1_updateflags(flagsp, (OPJ_UINT32)v, t1->flags_stride);
        }
        *flagsp |= T1_VISIT;
    }
}


void yl_opj_t1_enc_sigpass_step(   opj_t1_t *t1,
                                   opj_flag_t *flagsp,
                                   OPJ_INT32 *datap,
                                   OPJ_INT32 *datap_dec,
                                   OPJ_UINT32 orient,
                                   OPJ_INT32 bpno,
                                   OPJ_INT32 one,
                                   OPJ_INT32 *nmsedec,
                                   OPJ_BYTE type,
                                   OPJ_UINT32 vsc
                               )
{
    OPJ_INT32 v,v_dec;
    OPJ_UINT32 flag;

    opj_mqc_t *mqc = t1->mqc;	/* MQC component */

    flag = vsc ? (OPJ_UINT32)((*flagsp) & (~(T1_SIG_S | T1_SIG_SE | T1_SIG_SW | T1_SGN_S))) : (OPJ_UINT32)(*flagsp);
    if ((flag & T1_SIG_OTH) && !(flag & (T1_SIG | T1_VISIT)))
    {
        v = opj_int_abs(*datap) & one ? 1 : 0;
        if (v)
        {
            //datap_dec[0] =  (opj_int_abs(*datap)==(*datap) ? 1 : -1)* ( opj_int_abs( datap_dec[0] ) |  one);
            datap_dec[0] =  ( opj_int_abs( datap_dec[0] ) |  one);
            datap_dec[0] =  ( opj_int_abs( datap_dec[0] ) |   (one >> 1) );
        }
        else
        {
            v_dec = datap_dec[0] & one;
            datap_dec[0] =  ( opj_int_abs( datap_dec[0] ) &  (~one) );
            datap_dec[0] =  ( opj_int_abs( datap_dec[0] ) |   (one >> 1) );
        }

        datap_dec[0]=(*datap < 0 ? -1 : 1)*datap_dec[0]; //Sign bits coded when sample become SIGNIFICANT


        opj_mqc_setcurctx(mqc, opj_t1_getctxno_zc(flag, orient));	/* ESSAI */
        if (type == T1_TYPE_RAW)  	/* BYPASS/LAZY MODE */
        {
            opj_mqc_bypass_enc(mqc, (OPJ_UINT32)v);
        }
        else
        {
            opj_mqc_encode(mqc, (OPJ_UINT32)v);
        }
        if (v)
        {
            v = *datap < 0 ? 1 : 0;
            *nmsedec +=	opj_t1_getnmsedec_sig((OPJ_UINT32)opj_int_abs(*datap), (OPJ_UINT32)(bpno + T1_NMSEDEC_FRACBITS));
            opj_mqc_setcurctx(mqc, opj_t1_getctxno_sc(flag));	/* ESSAI */
            if (type == T1_TYPE_RAW)  	/* BYPASS/LAZY MODE */
            {
                opj_mqc_bypass_enc(mqc, (OPJ_UINT32)v);
            }
            else
            {
                opj_mqc_encode(mqc, (OPJ_UINT32)(v ^ opj_t1_getspb((OPJ_UINT32)flag)));//encode-sign
            }
            opj_t1_updateflags(flagsp, (OPJ_UINT32)v, t1->flags_stride);
        }
        *flagsp |= T1_VISIT;
    }
}



static INLINE void opj_t1_dec_sigpass_step_raw(
    opj_t1_t *t1,
    opj_flag_t *flagsp,
    OPJ_INT32 *datap,
    OPJ_INT32 orient,
    OPJ_INT32 oneplushalf,
    OPJ_INT32 vsc)
{
    OPJ_INT32 v, flag;
    opj_raw_t *raw = t1->raw;       /* RAW component */
    OPJ_ARG_NOT_USED(orient);

    flag = vsc ? ((*flagsp) & (~(T1_SIG_S | T1_SIG_SE | T1_SIG_SW | T1_SGN_S))) : (*flagsp);
    if ((flag & T1_SIG_OTH) && !(flag & (T1_SIG | T1_VISIT)))
    {
        if (opj_raw_decode(raw))
        {
            v = (OPJ_INT32)opj_raw_decode(raw);    /* ESSAI */
            *datap = v ? -oneplushalf : oneplushalf;
            opj_t1_updateflags(flagsp, (OPJ_UINT32)v, t1->flags_stride);
        }
        *flagsp |= T1_VISIT;
    }
}

INLINE void opj_t1_dec_sigpass_step_mqc(
    opj_t1_t *t1,
    opj_flag_t *flagsp,
    OPJ_INT32 *datap,
    OPJ_INT32 orient,
    OPJ_INT32 oneplushalf)
{
    OPJ_INT32 v, flag;

    opj_mqc_t *mqc = t1->mqc;       /* MQC component */

    flag = *flagsp;
    if ((flag & T1_SIG_OTH) && !(flag & (T1_SIG | T1_VISIT)))
    {
        opj_mqc_setcurctx(mqc, opj_t1_getctxno_zc((OPJ_UINT32)flag, (OPJ_UINT32)orient));
        if (opj_mqc_decode(mqc))
        {
            opj_mqc_setcurctx(mqc, opj_t1_getctxno_sc((OPJ_UINT32)flag));
            v = opj_mqc_decode(mqc) ^ opj_t1_getspb((OPJ_UINT32)flag);
            *datap = v ? -oneplushalf : oneplushalf;
            opj_t1_updateflags(flagsp, (OPJ_UINT32)v, t1->flags_stride);
        }
        *flagsp |= T1_VISIT;
    }
}                               /* VSC and  BYPASS by Antonin */

INLINE void opj_t1_dec_sigpass_step_mqc_vsc(
    opj_t1_t *t1,
    opj_flag_t *flagsp,
    OPJ_INT32 *datap,
    OPJ_INT32 orient,
    OPJ_INT32 oneplushalf,
    OPJ_INT32 vsc)
{
    OPJ_INT32 v, flag;

    opj_mqc_t *mqc = t1->mqc;       /* MQC component */

    flag = vsc ? ((*flagsp) & (~(T1_SIG_S | T1_SIG_SE | T1_SIG_SW | T1_SGN_S))) : (*flagsp);
    if ((flag & T1_SIG_OTH) && !(flag & (T1_SIG | T1_VISIT)))
    {
        opj_mqc_setcurctx(mqc, opj_t1_getctxno_zc((OPJ_UINT32)flag, (OPJ_UINT32)orient));
        if (opj_mqc_decode(mqc))
        {
            opj_mqc_setcurctx(mqc, opj_t1_getctxno_sc((OPJ_UINT32)flag));
            v = opj_mqc_decode(mqc) ^ opj_t1_getspb((OPJ_UINT32)flag);
            *datap = v ? -oneplushalf : oneplushalf;
            opj_t1_updateflags(flagsp, (OPJ_UINT32)v, t1->flags_stride);
        }
        *flagsp |= T1_VISIT;
    }
}                               /* VSC and  BYPASS by Antonin */



void opj_t1_enc_sigpass(opj_t1_t *t1,
                        OPJ_INT32 bpno,
                        OPJ_UINT32 orient,
                        OPJ_INT32 *nmsedec,
                        OPJ_BYTE type,
                        OPJ_UINT32 cblksty
                       )
{
    OPJ_UINT32 i, j, k, vsc;
    OPJ_INT32 one;

    *nmsedec = 0;
    one = 1 << (bpno + T1_NMSEDEC_FRACBITS);
    for (k = 0; k < t1->h; k += 4)
    {
        for (i = 0; i < t1->w; ++i)
        {
            for (j = k; j < k + 4 && j < t1->h; ++j)
            {
                vsc = ((cblksty & J2K_CCP_CBLKSTY_VSC) && (j == k + 3 || j == t1->h - 1)) ? 1 : 0;
                opj_t1_enc_sigpass_step(
                    t1,
                    &t1->flags[((j+1) * t1->flags_stride) + i + 1],
                    &t1->data[(j * t1->w) + i],
                    orient,
                    bpno,
                    one,
                    nmsedec,
                    type,
                    vsc);
            }
        }
    }
}

void yl_opj_t1_enc_sigpass(opj_t1_t *t1,
                           OPJ_INT32 * datap_dec,
                           OPJ_INT32 bpno,
                           OPJ_UINT32 orient,
                           OPJ_INT32 *nmsedec,
                           OPJ_BYTE type,
                           OPJ_UINT32 cblksty
                          )
{
    OPJ_UINT32 i, j, k, vsc;
    OPJ_INT32 one;

    *nmsedec = 0;
    one = 1 << (bpno + T1_NMSEDEC_FRACBITS);
    for (k = 0; k < t1->h; k += 4)
    {
        for (i = 0; i < t1->w; ++i)
        {
            for (j = k; j < k + 4 && j < t1->h; ++j)
            {
                vsc = ((cblksty & J2K_CCP_CBLKSTY_VSC) && (j == k + 3 || j == t1->h - 1)) ? 1 : 0;
                yl_opj_t1_enc_sigpass_step(
                    t1,
                    &t1->flags[((j+1) * t1->flags_stride) + i + 1],
                    &t1->data[(j * t1->w) + i],
                    &datap_dec[(j * t1->w) + i],
                    orient,
                    bpno,
                    one,
                    nmsedec,
                    type,
                    vsc);


            }
        }
    }
}


void opj_t1_dec_sigpass_raw(
    opj_t1_t *t1,
    OPJ_INT32 bpno,
    OPJ_INT32 orient,
    OPJ_INT32 cblksty)
{
    OPJ_INT32 one, half, oneplushalf, vsc;
    OPJ_UINT32 i, j, k;
    one = 1 << bpno;
    half = one >> 1;
    oneplushalf = one | half;
    for (k = 0; k < t1->h; k += 4)
    {
        for (i = 0; i < t1->w; ++i)
        {
            for (j = k; j < k + 4 && j < t1->h; ++j)
            {
                vsc = ((cblksty & J2K_CCP_CBLKSTY_VSC) && (j == k + 3 || j == t1->h - 1)) ? 1 : 0;
                opj_t1_dec_sigpass_step_raw(
                    t1,
                    &t1->flags[((j+1) * t1->flags_stride) + i + 1],
                    &t1->data[(j * t1->w) + i],
                    orient,
                    oneplushalf,
                    vsc);
            }
        }
    }
}                               /* VSC and  BYPASS by Antonin */

void opj_t1_dec_sigpass_mqc(
    opj_t1_t *t1,
    OPJ_INT32 bpno,
    OPJ_INT32 orient)
{
    OPJ_INT32 one, half, oneplushalf;
    OPJ_UINT32 i, j, k;
    OPJ_INT32 *data1 = t1->data;
    opj_flag_t *flags1 = &t1->flags[1];
    one = 1 << bpno;
    half = one >> 1;
    oneplushalf = one | half;
    for (k = 0; k < (t1->h & ~3u); k += 4)
    {
        for (i = 0; i < t1->w; ++i)
        {
            OPJ_INT32 *data2 = data1 + i;
            opj_flag_t *flags2 = flags1 + i;
            flags2 += t1->flags_stride;
            opj_t1_dec_sigpass_step_mqc(t1, flags2, data2, orient, oneplushalf);
            data2 += t1->w;
            flags2 += t1->flags_stride;
            opj_t1_dec_sigpass_step_mqc(t1, flags2, data2, orient, oneplushalf);
            data2 += t1->w;
            flags2 += t1->flags_stride;
            opj_t1_dec_sigpass_step_mqc(t1, flags2, data2, orient, oneplushalf);
            data2 += t1->w;
            flags2 += t1->flags_stride;
            opj_t1_dec_sigpass_step_mqc(t1, flags2, data2, orient, oneplushalf);
            data2 += t1->w;
        }
        data1 += t1->w << 2;
        flags1 += t1->flags_stride << 2;
    }
    for (i = 0; i < t1->w; ++i)
    {
        OPJ_INT32 *data2 = data1 + i;
        opj_flag_t *flags2 = flags1 + i;
        for (j = k; j < t1->h; ++j)
        {
            flags2 += t1->flags_stride;
            opj_t1_dec_sigpass_step_mqc(t1, flags2, data2, orient, oneplushalf);
            data2 += t1->w;
        }
    }
}                               /* VSC and  BYPASS by Antonin */

void opj_t1_dec_sigpass_mqc_vsc(
    opj_t1_t *t1,
    OPJ_INT32 bpno,
    OPJ_INT32 orient)
{
    OPJ_INT32 one, half, oneplushalf, vsc;
    OPJ_UINT32 i, j, k;
    one = 1 << bpno;
    half = one >> 1;
    oneplushalf = one | half;
    for (k = 0; k < t1->h; k += 4)
    {
        for (i = 0; i < t1->w; ++i)
        {
            for (j = k; j < k + 4 && j < t1->h; ++j)
            {
                vsc = (j == k + 3 || j == t1->h - 1) ? 1 : 0;
                opj_t1_dec_sigpass_step_mqc_vsc(
                    t1,
                    &t1->flags[((j+1) * t1->flags_stride) + i + 1],
                    &t1->data[(j * t1->w) + i],
                    orient,
                    oneplushalf,
                    vsc);
            }
        }
    }
}                               /* VSC and  BYPASS by Antonin */



void opj_t1_enc_refpass_step(   opj_t1_t *t1,
                                opj_flag_t *flagsp,
                                OPJ_INT32 *datap,
                                OPJ_INT32 bpno,
                                OPJ_INT32 one,
                                OPJ_INT32 *nmsedec,
                                OPJ_BYTE type,
                                OPJ_UINT32 vsc)
{
    OPJ_INT32 v;
    OPJ_UINT32 flag;

    opj_mqc_t *mqc = t1->mqc;	/* MQC component */

    flag = vsc ? (OPJ_UINT32)((*flagsp) & (~(T1_SIG_S | T1_SIG_SE | T1_SIG_SW | T1_SGN_S))) : (OPJ_UINT32)(*flagsp);
    if ((flag & (T1_SIG | T1_VISIT)) == T1_SIG)
    {
        *nmsedec += opj_t1_getnmsedec_ref((OPJ_UINT32)opj_int_abs(*datap), (OPJ_UINT32)(bpno + T1_NMSEDEC_FRACBITS));
        v = opj_int_abs(*datap) & one ? 1 : 0;
        opj_mqc_setcurctx(mqc, opj_t1_getctxno_mag(flag));	/* ESSAI */
        if (type == T1_TYPE_RAW)  	/* BYPASS/LAZY MODE */
        {
            opj_mqc_bypass_enc(mqc, (OPJ_UINT32)v);
        }
        else
        {
            opj_mqc_encode(mqc, (OPJ_UINT32)v);
        }
        *flagsp |= T1_REFINE;
    }
}



void yl_opj_t1_enc_refpass_step( opj_t1_t *t1,
                                 opj_flag_t *flagsp,
                                 OPJ_INT32 *datap,
                                 OPJ_INT32 *datap_dec,
                                 OPJ_INT32 bpno,
                                 OPJ_INT32 one,
                                 OPJ_INT32 *nmsedec,
                                 OPJ_BYTE type,
                                 OPJ_UINT32 vsc)
{
    OPJ_INT32 v,v_dec;
    OPJ_UINT32 flag;

    opj_mqc_t *mqc = t1->mqc;	/* MQC component */

    flag = vsc ? (OPJ_UINT32)((*flagsp) & (~(T1_SIG_S | T1_SIG_SE | T1_SIG_SW | T1_SGN_S))) : (OPJ_UINT32)(*flagsp);
    if ((flag & (T1_SIG | T1_VISIT)) == T1_SIG)
    {
        *nmsedec += opj_t1_getnmsedec_ref((OPJ_UINT32)opj_int_abs(*datap), (OPJ_UINT32)(bpno + T1_NMSEDEC_FRACBITS));
        v = opj_int_abs(*datap) & one ? 1 : 0;

        //Sign is coded in the other two passes befoer entering this passes, otherwise this point will not be reached
        if (v)
        {
            datap_dec[0] =  (opj_int_abs(*datap)==(*datap) ? 1 : -1)* ( opj_int_abs( datap_dec[0] ) |  one);
            datap_dec[0] =  (opj_int_abs(*datap)==(*datap) ? 1 : -1)* ( opj_int_abs( datap_dec[0] ) |  (one >> 1) );  //mid-point recon
        }
        else
        {
            v_dec=opj_int_abs( datap_dec[0] ) &  one;
            datap_dec[0] =  (opj_int_abs(*datap)==(*datap) ? 1 : -1)* ( opj_int_abs( datap_dec[0] ) &  (~ one) );
            datap_dec[0] =  (opj_int_abs(*datap)==(*datap) ? 1 : -1)* ( opj_int_abs( datap_dec[0] ) |  (one >> 1) );  //mid-point recon

        }

        opj_mqc_setcurctx(mqc, opj_t1_getctxno_mag(flag));	/* ESSAI */
        if (type == T1_TYPE_RAW)  	/* BYPASS/LAZY MODE */
        {
            opj_mqc_bypass_enc(mqc, (OPJ_UINT32)v);
        }
        else
        {
            opj_mqc_encode(mqc, (OPJ_UINT32)v);
        }
        *flagsp |= T1_REFINE;
    }
}



INLINE void opj_t1_dec_refpass_step_raw(
    opj_t1_t *t1,
    opj_flag_t *flagsp,
    OPJ_INT32 *datap,
    OPJ_INT32 poshalf,
    OPJ_INT32 neghalf,
    OPJ_INT32 vsc)
{
    OPJ_INT32 v, t, flag;

    opj_raw_t *raw = t1->raw;       /* RAW component */

    flag = vsc ? ((*flagsp) & (~(T1_SIG_S | T1_SIG_SE | T1_SIG_SW | T1_SGN_S))) : (*flagsp);
    if ((flag & (T1_SIG | T1_VISIT)) == T1_SIG)
    {
        v = (OPJ_INT32)opj_raw_decode(raw);
        t = v ? poshalf : neghalf;
        *datap += *datap < 0 ? -t : t;
        *flagsp |= T1_REFINE;
    }
}                               /* VSC and  BYPASS by Antonin  */

INLINE void opj_t1_dec_refpass_step_mqc(
    opj_t1_t *t1,
    opj_flag_t *flagsp,
    OPJ_INT32 *datap,
    OPJ_INT32 poshalf,
    OPJ_INT32 neghalf)
{
    OPJ_INT32 v, t, flag;

    opj_mqc_t *mqc = t1->mqc;       /* MQC component */

    flag = *flagsp;
    if ((flag & (T1_SIG | T1_VISIT)) == T1_SIG)
    {
        opj_mqc_setcurctx(mqc, opj_t1_getctxno_mag((OPJ_UINT32)flag));      /* ESSAI */
        v = opj_mqc_decode(mqc);
        t = v ? poshalf : neghalf;
        *datap += *datap < 0 ? -t : t;
        *flagsp |= T1_REFINE;
    }
}                               /* VSC and  BYPASS by Antonin  */

INLINE void opj_t1_dec_refpass_step_mqc_vsc(
    opj_t1_t *t1,
    opj_flag_t *flagsp,
    OPJ_INT32 *datap,
    OPJ_INT32 poshalf,
    OPJ_INT32 neghalf,
    OPJ_INT32 vsc)
{
    OPJ_INT32 v, t, flag;

    opj_mqc_t *mqc = t1->mqc;       /* MQC component */

    flag = vsc ? ((*flagsp) & (~(T1_SIG_S | T1_SIG_SE | T1_SIG_SW | T1_SGN_S))) : (*flagsp);
    if ((flag & (T1_SIG | T1_VISIT)) == T1_SIG)
    {
        opj_mqc_setcurctx(mqc, opj_t1_getctxno_mag((OPJ_UINT32)flag));      /* ESSAI */
        v = opj_mqc_decode(mqc);
        t = v ? poshalf : neghalf;
        *datap += *datap < 0 ? -t : t;
        *flagsp |= T1_REFINE;
    }
}                               /* VSC and  BYPASS by Antonin  */


void opj_t1_enc_refpass(
    opj_t1_t *t1,
    OPJ_INT32 bpno,
    OPJ_INT32 *nmsedec,
    OPJ_BYTE type,
    OPJ_UINT32 cblksty)
{
    OPJ_UINT32 i, j, k, vsc;
    OPJ_INT32 one;

    *nmsedec = 0;
    one = 1 << (bpno + T1_NMSEDEC_FRACBITS);
    for (k = 0; k < t1->h; k += 4)
    {
        for (i = 0; i < t1->w; ++i)
        {
            for (j = k; j < k + 4 && j < t1->h; ++j)
            {
                vsc = ((cblksty & J2K_CCP_CBLKSTY_VSC) && (j == k + 3 || j == t1->h - 1)) ? 1 : 0;
                opj_t1_enc_refpass_step(
                    t1,
                    &t1->flags[((j+1) * t1->flags_stride) + i + 1],
                    &t1->data[(j * t1->w) + i],
                    bpno,
                    one,
                    nmsedec,
                    type,
                    vsc);
            }
        }
    }
}


void yl_opj_t1_enc_refpass(
    opj_t1_t *t1,
    OPJ_INT32 *datap_dec,
    OPJ_INT32 bpno,
    OPJ_INT32 *nmsedec,
    OPJ_BYTE type,
    OPJ_UINT32 cblksty)
{
    OPJ_UINT32 i, j, k, vsc;
    OPJ_INT32 one;

    *nmsedec = 0;
    one = 1 << (bpno + T1_NMSEDEC_FRACBITS);
    for (k = 0; k < t1->h; k += 4)
    {
        for (i = 0; i < t1->w; ++i)
        {
            for (j = k; j < k + 4 && j < t1->h; ++j)
            {
                vsc = ((cblksty & J2K_CCP_CBLKSTY_VSC) && (j == k + 3 || j == t1->h - 1)) ? 1 : 0;
                yl_opj_t1_enc_refpass_step(
                    t1,
                    &t1->flags[((j+1) * t1->flags_stride) + i + 1],
                    &t1->data[(j * t1->w) + i],
                    &datap_dec[(j * t1->w) + i],
                    bpno,
                    one,
                    nmsedec,
                    type,
                    vsc);



            }
        }
    }
}


void opj_t1_dec_refpass_raw(
    opj_t1_t *t1,
    OPJ_INT32 bpno,
    OPJ_INT32 cblksty)
{
    OPJ_INT32 one, poshalf, neghalf;
    OPJ_UINT32 i, j, k;
    OPJ_INT32 vsc;
    one = 1 << bpno;
    poshalf = one >> 1;
    neghalf = bpno > 0 ? -poshalf : -1;
    for (k = 0; k < t1->h; k += 4)
    {
        for (i = 0; i < t1->w; ++i)
        {
            for (j = k; j < k + 4 && j < t1->h; ++j)
            {
                vsc = ((cblksty & J2K_CCP_CBLKSTY_VSC) && (j == k + 3 || j == t1->h - 1)) ? 1 : 0;
                opj_t1_dec_refpass_step_raw(
                    t1,
                    &t1->flags[((j+1) * t1->flags_stride) + i + 1],
                    &t1->data[(j * t1->w) + i],
                    poshalf,
                    neghalf,
                    vsc);
            }
        }
    }
}                               /* VSC and  BYPASS by Antonin */

void opj_t1_dec_refpass_mqc(
    opj_t1_t *t1,
    OPJ_INT32 bpno)
{
    OPJ_INT32 one, poshalf, neghalf;
    OPJ_UINT32 i, j, k;
    OPJ_INT32 *data1 = t1->data;
    opj_flag_t *flags1 = &t1->flags[1];
    one = 1 << bpno;
    poshalf = one >> 1;
    neghalf = bpno > 0 ? -poshalf : -1;
    for (k = 0; k < (t1->h & ~3u); k += 4)
    {
        for (i = 0; i < t1->w; ++i)
        {
            OPJ_INT32 *data2 = data1 + i;
            opj_flag_t *flags2 = flags1 + i;
            flags2 += t1->flags_stride;
            opj_t1_dec_refpass_step_mqc(t1, flags2, data2, poshalf, neghalf);
            data2 += t1->w;
            flags2 += t1->flags_stride;
            opj_t1_dec_refpass_step_mqc(t1, flags2, data2, poshalf, neghalf);
            data2 += t1->w;
            flags2 += t1->flags_stride;
            opj_t1_dec_refpass_step_mqc(t1, flags2, data2, poshalf, neghalf);
            data2 += t1->w;
            flags2 += t1->flags_stride;
            opj_t1_dec_refpass_step_mqc(t1, flags2, data2, poshalf, neghalf);
            data2 += t1->w;
        }
        data1 += t1->w << 2;
        flags1 += t1->flags_stride << 2;
    }
    for (i = 0; i < t1->w; ++i)
    {
        OPJ_INT32 *data2 = data1 + i;
        opj_flag_t *flags2 = flags1 + i;
        for (j = k; j < t1->h; ++j)
        {
            flags2 += t1->flags_stride;
            opj_t1_dec_refpass_step_mqc(t1, flags2, data2, poshalf, neghalf);
            data2 += t1->w;
        }
    }
}                               /* VSC and  BYPASS by Antonin */

void opj_t1_dec_refpass_mqc_vsc(
    opj_t1_t *t1,
    OPJ_INT32 bpno)
{
    OPJ_INT32 one, poshalf, neghalf;
    OPJ_UINT32 i, j, k;
    OPJ_INT32 vsc;
    one = 1 << bpno;
    poshalf = one >> 1;
    neghalf = bpno > 0 ? -poshalf : -1;
    for (k = 0; k < t1->h; k += 4)
    {
        for (i = 0; i < t1->w; ++i)
        {
            for (j = k; j < k + 4 && j < t1->h; ++j)
            {
                vsc = ((j == k + 3 || j == t1->h - 1)) ? 1 : 0;
                opj_t1_dec_refpass_step_mqc_vsc(
                    t1,
                    &t1->flags[((j+1) * t1->flags_stride) + i + 1],
                    &t1->data[(j * t1->w) + i],
                    poshalf,
                    neghalf,
                    vsc);
            }
        }
    }
}                               /* VSC and  BYPASS by Antonin */


void opj_t1_enc_clnpass_step(
    opj_t1_t *t1,
    opj_flag_t *flagsp,
    OPJ_INT32 *datap,
    OPJ_UINT32 orient,
    OPJ_INT32 bpno,
    OPJ_INT32 one,
    OPJ_INT32 *nmsedec,
    OPJ_UINT32 partial,
    OPJ_UINT32 vsc)
{
    OPJ_INT32 v;
    OPJ_UINT32 flag;

    opj_mqc_t *mqc = t1->mqc;	/* MQC component */

    flag = vsc ? (OPJ_UINT32)((*flagsp) & (~(T1_SIG_S | T1_SIG_SE | T1_SIG_SW | T1_SGN_S))) : (OPJ_UINT32)(*flagsp);
    if (partial)
    {
        goto LABEL_PARTIAL;
    }
    if (!(*flagsp & (T1_SIG | T1_VISIT)))
    {
        opj_mqc_setcurctx(mqc, opj_t1_getctxno_zc(flag, orient));
        v = opj_int_abs(*datap) & one ? 1 : 0;
        opj_mqc_encode(mqc, (OPJ_UINT32)v);
        if (v)
        {
LABEL_PARTIAL:
            *nmsedec += opj_t1_getnmsedec_sig((OPJ_UINT32)opj_int_abs(*datap), (OPJ_UINT32)(bpno + T1_NMSEDEC_FRACBITS));
            opj_mqc_setcurctx(mqc, opj_t1_getctxno_sc(flag));
            v = *datap < 0 ? 1 : 0;
            opj_mqc_encode(mqc, (OPJ_UINT32)(v ^ opj_t1_getspb((OPJ_UINT32)flag)));
            opj_t1_updateflags(flagsp, (OPJ_UINT32)v, t1->flags_stride);
        }
    }
    *flagsp &= ~T1_VISIT;
}

static void opj_t1_dec_clnpass_step_partial(
    opj_t1_t *t1,
    opj_flag_t *flagsp,
    OPJ_INT32 *datap,
    OPJ_INT32 orient,
    OPJ_INT32 oneplushalf)
{
    OPJ_INT32 v, flag;
    opj_mqc_t *mqc = t1->mqc;	/* MQC component */

    OPJ_ARG_NOT_USED(orient);

    flag = *flagsp;
    opj_mqc_setcurctx(mqc, opj_t1_getctxno_sc((OPJ_UINT32)flag));
    v = opj_mqc_decode(mqc) ^ opj_t1_getspb((OPJ_UINT32)flag);
    *datap = v ? -oneplushalf : oneplushalf;
    opj_t1_updateflags(flagsp, (OPJ_UINT32)v, t1->flags_stride);
    *flagsp &= ~T1_VISIT;
}				/* VSC and  BYPASS by Antonin */

static void opj_t1_dec_clnpass_step(
    opj_t1_t *t1,
    opj_flag_t *flagsp,
    OPJ_INT32 *datap,
    OPJ_INT32 orient,
    OPJ_INT32 oneplushalf)
{
    OPJ_INT32 v, flag;

    opj_mqc_t *mqc = t1->mqc;	/* MQC component */

    flag = *flagsp;
    if (!(flag & (T1_SIG | T1_VISIT)))
    {
        opj_mqc_setcurctx(mqc, opj_t1_getctxno_zc((OPJ_UINT32)flag, (OPJ_UINT32)orient));
        if (opj_mqc_decode(mqc))
        {
            opj_mqc_setcurctx(mqc, opj_t1_getctxno_sc((OPJ_UINT32)flag));
            v = opj_mqc_decode(mqc) ^ opj_t1_getspb((OPJ_UINT32)flag);
            *datap = v ? -oneplushalf : oneplushalf;
            opj_t1_updateflags(flagsp, (OPJ_UINT32)v, t1->flags_stride);
        }
    }
    *flagsp &= ~T1_VISIT;
}				/* VSC and  BYPASS by Antonin */

static void opj_t1_dec_clnpass_step_vsc(
    opj_t1_t *t1,
    opj_flag_t *flagsp,
    OPJ_INT32 *datap,
    OPJ_INT32 orient,
    OPJ_INT32 oneplushalf,
    OPJ_INT32 partial,
    OPJ_INT32 vsc)
{
    OPJ_INT32 v, flag;

    opj_mqc_t *mqc = t1->mqc;	/* MQC component */

    flag = vsc ? ((*flagsp) & (~(T1_SIG_S | T1_SIG_SE | T1_SIG_SW | T1_SGN_S))) : (*flagsp);
    if (partial)
    {
        goto LABEL_PARTIAL;
    }
    if (!(flag & (T1_SIG | T1_VISIT)))
    {
        opj_mqc_setcurctx(mqc, opj_t1_getctxno_zc((OPJ_UINT32)flag, (OPJ_UINT32)orient));
        if (opj_mqc_decode(mqc))
        {
LABEL_PARTIAL:
            opj_mqc_setcurctx(mqc, opj_t1_getctxno_sc((OPJ_UINT32)flag));
            v = opj_mqc_decode(mqc) ^ opj_t1_getspb((OPJ_UINT32)flag);
            *datap = v ? -oneplushalf : oneplushalf;
            opj_t1_updateflags(flagsp, (OPJ_UINT32)v, t1->flags_stride);
        }
    }
    *flagsp &= ~T1_VISIT;
}

void opj_t1_enc_clnpass(
    opj_t1_t *t1,
    OPJ_INT32 bpno,
    OPJ_UINT32 orient,
    OPJ_INT32 *nmsedec,
    OPJ_UINT32 cblksty)
{
    OPJ_UINT32 i, j, k;
    OPJ_INT32 one;
    OPJ_UINT32 agg, runlen, vsc;

    opj_mqc_t *mqc = t1->mqc;	/* MQC component */

    *nmsedec = 0;
    one = 1 << (bpno + T1_NMSEDEC_FRACBITS);
    for (k = 0; k < t1->h; k += 4)
    {
        for (i = 0; i < t1->w; ++i)
        {
            if (k + 3 < t1->h)
            {
                if (cblksty & J2K_CCP_CBLKSTY_VSC)
                {
                    agg = !(MACRO_t1_flags(1 + k,1 + i) & (T1_SIG | T1_VISIT | T1_SIG_OTH)
                            || MACRO_t1_flags(1 + k + 1,1 + i) & (T1_SIG | T1_VISIT | T1_SIG_OTH)
                            || MACRO_t1_flags(1 + k + 2,1 + i) & (T1_SIG | T1_VISIT | T1_SIG_OTH)
                            || (MACRO_t1_flags(1 + k + 3,1 + i)
                                & (~(T1_SIG_S | T1_SIG_SE | T1_SIG_SW |	T1_SGN_S))) & (T1_SIG | T1_VISIT | T1_SIG_OTH));
                }
                else
                {
                    agg = !(MACRO_t1_flags(1 + k,1 + i) & (T1_SIG | T1_VISIT | T1_SIG_OTH)
                            || MACRO_t1_flags(1 + k + 1,1 + i) & (T1_SIG | T1_VISIT | T1_SIG_OTH)
                            || MACRO_t1_flags(1 + k + 2,1 + i) & (T1_SIG | T1_VISIT | T1_SIG_OTH)
                            || MACRO_t1_flags(1 + k + 3,1 + i) & (T1_SIG | T1_VISIT | T1_SIG_OTH));
                }
            }
            else
            {
                agg = 0;
            }
            if (agg)
            {
                for (runlen = 0; runlen < 4; ++runlen)
                {
                    if (opj_int_abs(t1->data[((k + runlen)*t1->w) + i]) & one)
                        break;
                }
                opj_mqc_setcurctx(mqc, T1_CTXNO_AGG);
                opj_mqc_encode(mqc, runlen != 4);
                if (runlen == 4)
                {
                    continue;
                }
                opj_mqc_setcurctx(mqc, T1_CTXNO_UNI);
                opj_mqc_encode(mqc, runlen >> 1);
                opj_mqc_encode(mqc, runlen & 1);
            }
            else
            {
                runlen = 0;
            }
            for (j = k + runlen; j < k + 4 && j < t1->h; ++j)
            {
                vsc = ((cblksty & J2K_CCP_CBLKSTY_VSC) && (j == k + 3 || j == t1->h - 1)) ? 1 : 0;
                opj_t1_enc_clnpass_step(
                    t1,
                    &t1->flags[((j+1) * t1->flags_stride) + i + 1],
                    &t1->data[(j * t1->w) + i],
                    orient,
                    bpno,
                    one,
                    nmsedec,
                    agg && (j == k + runlen),
                    vsc);
            }
        }
    }
}



void yl_opj_t1_enc_clnpass(
    opj_t1_t *t1,
    OPJ_INT32 *datap_dec,
    OPJ_INT32 bpno,
    OPJ_UINT32 orient,
    OPJ_INT32 *nmsedec,
    OPJ_UINT32 cblksty)
{
    OPJ_UINT32 i, j, k;
    OPJ_INT32 one;
    OPJ_UINT32 agg, runlen, vsc;

    opj_mqc_t *mqc = t1->mqc;	/* MQC component */

    *nmsedec = 0;
    one = 1 << (bpno + T1_NMSEDEC_FRACBITS);
    for (k = 0; k < t1->h; k += 4)
    {
        for (i = 0; i < t1->w; ++i)
        {
            if (k + 3 < t1->h)
            {
                if (cblksty & J2K_CCP_CBLKSTY_VSC)
                {
                    agg = !(MACRO_t1_flags(1 + k,1 + i) & (T1_SIG | T1_VISIT | T1_SIG_OTH)
                            || MACRO_t1_flags(1 + k + 1,1 + i) & (T1_SIG | T1_VISIT | T1_SIG_OTH)
                            || MACRO_t1_flags(1 + k + 2,1 + i) & (T1_SIG | T1_VISIT | T1_SIG_OTH)
                            || (MACRO_t1_flags(1 + k + 3,1 + i)
                                & (~(T1_SIG_S | T1_SIG_SE | T1_SIG_SW |	T1_SGN_S))) & (T1_SIG | T1_VISIT | T1_SIG_OTH));
                }
                else
                {
                    agg = !(MACRO_t1_flags(1 + k,1 + i) & (T1_SIG | T1_VISIT | T1_SIG_OTH)
                            || MACRO_t1_flags(1 + k + 1,1 + i) & (T1_SIG | T1_VISIT | T1_SIG_OTH)
                            || MACRO_t1_flags(1 + k + 2,1 + i) & (T1_SIG | T1_VISIT | T1_SIG_OTH)
                            || MACRO_t1_flags(1 + k + 3,1 + i) & (T1_SIG | T1_VISIT | T1_SIG_OTH));
                }
            }
            else
            {
                agg = 0;
            }
            if (agg)
            {
                for (runlen = 0; runlen < 4; ++runlen)
                {
                    if (opj_int_abs(t1->data[((k + runlen)*t1->w) + i]) & one)
                        break;
                }
                opj_mqc_setcurctx(mqc, T1_CTXNO_AGG);
                opj_mqc_encode(mqc, runlen != 4);
                if (runlen == 4)
                {
                    continue;
                }
                opj_mqc_setcurctx(mqc, T1_CTXNO_UNI);
                opj_mqc_encode(mqc, runlen >> 1);
                opj_mqc_encode(mqc, runlen & 1);
            }
            else
            {
                runlen = 0;
            }
            for (j = k + runlen; j < k + 4 && j < t1->h; ++j)
            {
                vsc = ((cblksty & J2K_CCP_CBLKSTY_VSC) && (j == k + 3 || j == t1->h - 1)) ? 1 : 0;


                opj_t1_enc_clnpass_step(
                    t1,
                    &t1->flags[((j+1) * t1->flags_stride) + i + 1],
                    &t1->data[(j * t1->w) + i],
                    orient,
                    bpno,
                    one,
                    nmsedec,
                    agg && (j == k + runlen),
                    vsc);
            }
        }
    }


    //yl
    OPJ_INT32 v_dec;
    for (i = 0; i < t1->w; ++i)
    {
        for (j = 0; j < t1->h; ++j)
        {
            OPJ_INT32 v = opj_int_abs(t1->data[(j * t1->w) + i]) & one ? 1 : 0;
            if (v)
            {
                datap_dec[(j * t1->w) + i] = ( ( t1->data[(j * t1->w) + i] >0 ) ? 1 : -1)* ( opj_int_abs( datap_dec[(j * t1->w) + i] ) |  one);
                datap_dec[(j * t1->w) + i] = ( ( t1->data[(j * t1->w) + i] >0 ) ? 1 : -1)* ( opj_int_abs( datap_dec[(j * t1->w) + i] ) |  (one >> 1)   );

            }
            else
            {
                v_dec = datap_dec[(j * t1->w) + i] & one ;

            }
        }
    }
    //yl:end


}



static void opj_t1_dec_clnpass(
    opj_t1_t *t1,
    OPJ_INT32 bpno,
    OPJ_INT32 orient,
    OPJ_INT32 cblksty)
{
    OPJ_INT32 one, half, oneplushalf, agg, runlen, vsc;
    OPJ_UINT32 i, j, k;
    OPJ_INT32 segsym = cblksty & J2K_CCP_CBLKSTY_SEGSYM;

    opj_mqc_t *mqc = t1->mqc;	/* MQC component */

    one = 1 << bpno;
    half = one >> 1;
    oneplushalf = one | half;
    if (cblksty & J2K_CCP_CBLKSTY_VSC)
    {
        for (k = 0; k < t1->h; k += 4)
        {
            for (i = 0; i < t1->w; ++i)
            {
                if (k + 3 < t1->h)
                {
                    agg = !(MACRO_t1_flags(1 + k,1 + i) & (T1_SIG | T1_VISIT | T1_SIG_OTH)
                            || MACRO_t1_flags(1 + k + 1,1 + i) & (T1_SIG | T1_VISIT | T1_SIG_OTH)
                            || MACRO_t1_flags(1 + k + 2,1 + i) & (T1_SIG | T1_VISIT | T1_SIG_OTH)
                            || (MACRO_t1_flags(1 + k + 3,1 + i)
                                & (~(T1_SIG_S | T1_SIG_SE | T1_SIG_SW |	T1_SGN_S))) & (T1_SIG | T1_VISIT | T1_SIG_OTH));
                }
                else
                {
                    agg = 0;
                }
                if (agg)
                {
                    opj_mqc_setcurctx(mqc, T1_CTXNO_AGG);
                    if (!opj_mqc_decode(mqc))
                    {
                        continue;
                    }
                    opj_mqc_setcurctx(mqc, T1_CTXNO_UNI);
                    runlen = opj_mqc_decode(mqc);
                    runlen = (runlen << 1) | opj_mqc_decode(mqc);
                }
                else
                {
                    runlen = 0;
                }
                for (j = k + (OPJ_UINT32)runlen; j < k + 4 && j < t1->h; ++j)
                {
                    vsc = (j == k + 3 || j == t1->h - 1) ? 1 : 0;
                    opj_t1_dec_clnpass_step_vsc(
                        t1,
                        &t1->flags[((j+1) * t1->flags_stride) + i + 1],
                        &t1->data[(j * t1->w) + i],
                        orient,
                        oneplushalf,
                        agg && (j == k + (OPJ_UINT32)runlen),
                        vsc);
                }
            }
        }
    }
    else
    {
        OPJ_INT32 *data1 = t1->data;
        opj_flag_t *flags1 = &t1->flags[1];
        for (k = 0; k < (t1->h & ~3u); k += 4)
        {
            for (i = 0; i < t1->w; ++i)
            {
                OPJ_INT32 *data2 = data1 + i;
                opj_flag_t *flags2 = flags1 + i;
                agg = !(MACRO_t1_flags(1 + k,1 + i) & (T1_SIG | T1_VISIT | T1_SIG_OTH)
                        || MACRO_t1_flags(1 + k + 1,1 + i) & (T1_SIG | T1_VISIT | T1_SIG_OTH)
                        || MACRO_t1_flags(1 + k + 2,1 + i) & (T1_SIG | T1_VISIT | T1_SIG_OTH)
                        || MACRO_t1_flags(1 + k + 3,1 + i) & (T1_SIG | T1_VISIT | T1_SIG_OTH));
                if (agg)
                {
                    opj_mqc_setcurctx(mqc, T1_CTXNO_AGG);
                    if (!opj_mqc_decode(mqc))
                    {
                        continue;
                    }
                    opj_mqc_setcurctx(mqc, T1_CTXNO_UNI);
                    runlen = opj_mqc_decode(mqc);
                    runlen = (runlen << 1) | opj_mqc_decode(mqc);
                    flags2 += (OPJ_UINT32)runlen * t1->flags_stride;
                    data2 += (OPJ_UINT32)runlen * t1->w;
                    for (j = k + (OPJ_UINT32)runlen; j < k + 4 && j < t1->h; ++j)
                    {
                        flags2 += t1->flags_stride;
                        if (agg && (j == k + (OPJ_UINT32)runlen))
                        {
                            opj_t1_dec_clnpass_step_partial(t1, flags2, data2, orient, oneplushalf);
                        }
                        else
                        {
                            opj_t1_dec_clnpass_step(t1, flags2, data2, orient, oneplushalf);
                        }
                        data2 += t1->w;
                    }
                }
                else
                {
                    flags2 += t1->flags_stride;
                    opj_t1_dec_clnpass_step(t1, flags2, data2, orient, oneplushalf);
                    data2 += t1->w;
                    flags2 += t1->flags_stride;
                    opj_t1_dec_clnpass_step(t1, flags2, data2, orient, oneplushalf);
                    data2 += t1->w;
                    flags2 += t1->flags_stride;
                    opj_t1_dec_clnpass_step(t1, flags2, data2, orient, oneplushalf);
                    data2 += t1->w;
                    flags2 += t1->flags_stride;
                    opj_t1_dec_clnpass_step(t1, flags2, data2, orient, oneplushalf);
                    data2 += t1->w;
                }
            }
            data1 += t1->w << 2;
            flags1 += t1->flags_stride << 2;
        }
        for (i = 0; i < t1->w; ++i)
        {
            OPJ_INT32 *data2 = data1 + i;
            opj_flag_t *flags2 = flags1 + i;
            for (j = k; j < t1->h; ++j)
            {
                flags2 += t1->flags_stride;
                opj_t1_dec_clnpass_step(t1, flags2, data2, orient, oneplushalf);
                data2 += t1->w;
            }
        }
    }

    if (segsym)
    {
        OPJ_INT32 v = 0;
        opj_mqc_setcurctx(mqc, T1_CTXNO_UNI);
        v = opj_mqc_decode(mqc);
        v = (v << 1) | opj_mqc_decode(mqc);
        v = (v << 1) | opj_mqc_decode(mqc);
        v = (v << 1) | opj_mqc_decode(mqc);
        /*
        if (v!=0xa) {
        	opj_event_msg(t1->cinfo, EVT_WARNING, "Bad segmentation symbol %x\n", v);
        }
        */
    }
}				/* VSC and  BYPASS by Antonin */


/** mod fixed_quality */
static OPJ_FLOAT64 opj_t1_getwmsedec(
    OPJ_INT32 nmsedec,
    OPJ_UINT32 compno,
    OPJ_UINT32 level,
    OPJ_UINT32 orient,
    OPJ_INT32 bpno,
    OPJ_UINT32 qmfbid,
    OPJ_FLOAT64 stepsize,
    OPJ_UINT32 numcomps,
    const OPJ_FLOAT64 * mct_norms)
{
    OPJ_FLOAT64 w1 = 1, w2, wmsedec;
    OPJ_ARG_NOT_USED(numcomps);

    if (mct_norms)
    {
        w1 = mct_norms[compno];
    }

    if (qmfbid == 1)
    {
        w2 = opj_dwt_getnorm(level, orient);
    }
    else  	/* if (qmfbid == 0) */
    {
        w2 = opj_dwt_getnorm_real(level, orient);
    }

    wmsedec = w1 * w2 * stepsize * (1 << bpno);
    wmsedec *= wmsedec * nmsedec / 8192.0;

    return wmsedec;
}

OPJ_BOOL opj_t1_allocate_buffers(
    opj_t1_t *t1,
    OPJ_UINT32 w,
    OPJ_UINT32 h)
{
    OPJ_UINT32 datasize=w * h;
    OPJ_UINT32 flagssize;

    if(datasize > t1->datasize)
    {
        opj_aligned_free(t1->data);
        t1->data = (OPJ_INT32*) opj_aligned_malloc(datasize * sizeof(OPJ_INT32));
        if(!t1->data)
        {
            return OPJ_FALSE;
        }
        t1->datasize=datasize;
    }
    memset(t1->data,0,datasize * sizeof(OPJ_INT32));

    t1->flags_stride=w+2;
    flagssize=t1->flags_stride * (h+2);

    if(flagssize > t1->flagssize)
    {
        opj_aligned_free(t1->flags);
        t1->flags = (opj_flag_t*) opj_aligned_malloc(flagssize * sizeof(opj_flag_t));
        if(!t1->flags)
        {
            return OPJ_FALSE;
        }
        t1->flagssize=flagssize;
    }
    memset(t1->flags,0,flagssize * sizeof(opj_flag_t));

    t1->w=w;
    t1->h=h;

    return OPJ_TRUE;
}

/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */
/**
 * Creates a new Tier 1 handle
 * and initializes the look-up tables of the Tier-1 coder/decoder
 * @return a new T1 handle if successful, returns NULL otherwise
*/
opj_t1_t* opj_t1_create()
{
    opj_t1_t *l_t1 = 00;

    l_t1 = (opj_t1_t*) opj_malloc(sizeof(opj_t1_t));
    if (!l_t1)
    {
        return 00;
    }
    memset(l_t1,0,sizeof(opj_t1_t));

    /* create MQC and RAW handles */
    l_t1->mqc = opj_mqc_create();
    if (! l_t1->mqc)
    {
        opj_t1_destroy(l_t1);
        return 00;
    }

    l_t1->raw = opj_raw_create();
    if (! l_t1->raw)
    {
        opj_t1_destroy(l_t1);
        return 00;
    }

    return l_t1;
}


/**
 * Destroys a previously created T1 handle
 *
 * @param p_t1 Tier 1 handle to destroy
*/
void opj_t1_destroy(opj_t1_t *p_t1)
{
    if (! p_t1)
    {
        return;
    }

    /* destroy MQC and RAW handles */
    opj_mqc_destroy(p_t1->mqc);
    p_t1->mqc = 00;
    opj_raw_destroy(p_t1->raw);
    p_t1->raw = 00;

    if (p_t1->data)
    {
        opj_aligned_free(p_t1->data);
        p_t1->data = 00;
    }

    if (p_t1->flags)
    {
        opj_aligned_free(p_t1->flags);
        p_t1->flags = 00;
    }

    opj_free(p_t1);
}




OPJ_BOOL opj_t1_decode_cblks(   opj_t1_t* t1,
                                opj_tcd_tilecomp_t* tilec,
                                opj_tccp_t* tccp
                            )
{
    OPJ_UINT32 resno, bandno, precno, cblkno;
    OPJ_UINT32 tile_w = (OPJ_UINT32)(tilec->x1 - tilec->x0);

    for (resno = 0; resno < tilec->minimum_num_resolutions; ++resno)
    {
        opj_tcd_resolution_t* res = &tilec->resolutions[resno];

        for (bandno = 0; bandno < res->numbands; ++bandno)
        {
            opj_tcd_band_t* restrict band = &res->bands[bandno];

            for (precno = 0; precno < res->pw * res->ph; ++precno)
            {
                opj_tcd_precinct_t* precinct = &band->precincts[precno];

                for (cblkno = 0; cblkno < precinct->cw * precinct->ch; ++cblkno)
                {
                    opj_tcd_cblk_dec_t* cblk = &precinct->cblks.dec[cblkno];
                    OPJ_INT32* restrict datap;
                    /*void* restrict tiledp;*/
                    OPJ_UINT32 cblk_w, cblk_h;
                    OPJ_INT32 x, y;
                    OPJ_UINT32 i, j;

                    if (OPJ_FALSE == opj_t1_decode_cblk(
                                t1,
                                cblk,
                                band->bandno,
                                (OPJ_UINT32)tccp->roishift,
                                tccp->cblksty))
                    {
                        return OPJ_FALSE;
                    }

                    x = cblk->x0 - band->x0;
                    y = cblk->y0 - band->y0;
                    if (band->bandno & 1)
                    {
                        opj_tcd_resolution_t* pres = &tilec->resolutions[resno - 1];
                        x += pres->x1 - pres->x0;
                    }
                    if (band->bandno & 2)
                    {
                        opj_tcd_resolution_t* pres = &tilec->resolutions[resno - 1];
                        y += pres->y1 - pres->y0;
                    }

                    datap=t1->data;
                    cblk_w = t1->w;
                    cblk_h = t1->h;

                    if (tccp->roishift)
                    {
                        OPJ_INT32 thresh = 1 << tccp->roishift;
                        for (j = 0; j < cblk_h; ++j)
                        {
                            for (i = 0; i < cblk_w; ++i)
                            {
                                OPJ_INT32 val = datap[(j * cblk_w) + i];
                                OPJ_INT32 mag = abs(val);
                                if (mag >= thresh)
                                {
                                    mag >>= tccp->roishift;
                                    datap[(j * cblk_w) + i] = val < 0 ? -mag : mag;
                                }
                            }
                        }
                    }

                    /*tiledp=(void*)&tilec->data[(y * tile_w) + x];*/
                    if (tccp->qmfbid == 1)
                    {
                        OPJ_INT32* restrict tiledp = &tilec->data[(OPJ_UINT32)y * tile_w + (OPJ_UINT32)x];
                        for (j = 0; j < cblk_h; ++j)
                        {
                            for (i = 0; i < cblk_w; ++i)
                            {
                                OPJ_INT32 tmp = datap[(j * cblk_w) + i];
                                ((OPJ_INT32*)tiledp)[(j * tile_w) + i] = tmp / 2;
                            }
                        }
                    }
                    else  		/* if (tccp->qmfbid == 0) */
                    {
                        OPJ_FLOAT32* restrict tiledp = (OPJ_FLOAT32*) &tilec->data[(OPJ_UINT32)y * tile_w + (OPJ_UINT32)x];
                        for (j = 0; j < cblk_h; ++j)
                        {
                            OPJ_FLOAT32* restrict tiledp2 = tiledp;
                            for (i = 0; i < cblk_w; ++i)
                            {
                                OPJ_FLOAT32 tmp = (OPJ_FLOAT32)*datap * band->stepsize;

                                *tiledp2 = tmp;
                                datap++;
                                tiledp2++;

                            }
                            tiledp += tile_w;
                        }

                        //classfication


                    }
                    /*opj_free(cblk->data);
                    opj_free(cblk->segs);*/
                    /*cblk->segs = 00;*/
                } /* cblkno */
                /*opj_free(precinct->cblks.dec);*/
            } /* precno */
        } /* bandno */
    } /* resno */
    return OPJ_TRUE;
}


OPJ_BOOL opj_t1_decode_cblk(opj_t1_t *t1,
                            opj_tcd_cblk_dec_t* cblk,
                            OPJ_UINT32 orient,
                            OPJ_UINT32 roishift,
                            OPJ_UINT32 cblksty)
{
    opj_raw_t *raw = t1->raw;	/* RAW component */
    opj_mqc_t *mqc = t1->mqc;	/* MQC component */

    OPJ_INT32 bpno;
    OPJ_UINT32 passtype;
    OPJ_UINT32 segno, passno;
    OPJ_BYTE type = T1_TYPE_MQ; /* BYPASS mode */

    if(!opj_t1_allocate_buffers(
                t1,
                (OPJ_UINT32)(cblk->x1 - cblk->x0),
                (OPJ_UINT32)(cblk->y1 - cblk->y0)))
    {
        return OPJ_FALSE;
    }

    bpno = (OPJ_INT32)(roishift + cblk->numbps - 1);
    passtype = 2;

    opj_mqc_resetstates(mqc);
    opj_mqc_setstate(mqc, T1_CTXNO_UNI, 0, 46);
    opj_mqc_setstate(mqc, T1_CTXNO_AGG, 0, 3);
    opj_mqc_setstate(mqc, T1_CTXNO_ZC, 0, 4);

    for (segno = 0; segno < cblk->real_num_segs; ++segno)
    {
        opj_tcd_seg_t *seg = &cblk->segs[segno];

        /* BYPASS mode */
        type = ((bpno <= ((OPJ_INT32) (cblk->numbps) - 1) - 4) && (passtype < 2) && (cblksty & J2K_CCP_CBLKSTY_LAZY)) ? T1_TYPE_RAW : T1_TYPE_MQ;
        /* FIXME: slviewer gets here with a null pointer. Why? Partially downloaded and/or corrupt textures? */
        if(seg->data == 00)
        {
            continue;
        }
        if (type == T1_TYPE_RAW)
        {
            opj_raw_init_dec(raw, (*seg->data) + seg->dataindex, seg->len);
        }
        else
        {
            if (OPJ_FALSE == opj_mqc_init_dec(mqc, (*seg->data) + seg->dataindex, seg->len))
            {
                return OPJ_FALSE;
            }
        }

        for (passno = 0; passno < seg->real_num_passes; ++passno)
        {
            switch (passtype)
            {
            case 0:
                if (type == T1_TYPE_RAW)
                {
                    opj_t1_dec_sigpass_raw(t1, bpno+1, (OPJ_INT32)orient, (OPJ_INT32)cblksty);
                }
                else
                {
                    if (cblksty & J2K_CCP_CBLKSTY_VSC)
                    {
                        opj_t1_dec_sigpass_mqc_vsc(t1, bpno+1, (OPJ_INT32)orient);
                    }
                    else
                    {
                        opj_t1_dec_sigpass_mqc(t1, bpno+1, (OPJ_INT32)orient);
                    }
                }
                break;
            case 1:
                if (type == T1_TYPE_RAW)
                {
                    opj_t1_dec_refpass_raw(t1, bpno+1, (OPJ_INT32)cblksty);
                }
                else
                {
                    if (cblksty & J2K_CCP_CBLKSTY_VSC)
                    {
                        opj_t1_dec_refpass_mqc_vsc(t1, bpno+1);
                    }
                    else
                    {
                        opj_t1_dec_refpass_mqc(t1, bpno+1);
                    }
                }
                break;
            case 2:
                opj_t1_dec_clnpass(t1, bpno+1, (OPJ_INT32)orient, (OPJ_INT32)cblksty);
                break;
            }

            if ((cblksty & J2K_CCP_CBLKSTY_RESET) && type == T1_TYPE_MQ)
            {
                opj_mqc_resetstates(mqc);
                opj_mqc_setstate(mqc, T1_CTXNO_UNI, 0, 46);
                opj_mqc_setstate(mqc, T1_CTXNO_AGG, 0, 3);
                opj_mqc_setstate(mqc, T1_CTXNO_ZC, 0, 4);
            }
            if (++passtype == 3)
            {
                passtype = 0;
                bpno--;
            }
        }
    }
    return OPJ_TRUE;
}




OPJ_BOOL opj_t1_encode_cblks(   opj_t1_t *t1,
                                opj_tcd_tile_t *tile,
                                opj_tcp_t *tcp,
                                const OPJ_FLOAT64 * mct_norms
                            )
{
#ifdef IS_DEBUG_PRINTF
    FILE* fp_dwt_idx = fopen("dwt_idx.bin","wb+");
#endif
    OPJ_UINT32 compno, resno, bandno, precno, cblkno;

    tile->distotile = 0;		/* fixed_quality */

    for (compno = 0; compno < tile->numcomps; ++compno)
    {
        opj_tcd_tilecomp_t* tilec = &tile->comps[compno];
        opj_tccp_t* tccp = &tcp->tccps[compno];
        OPJ_UINT32 tile_w = (OPJ_UINT32)(tilec->x1 - tilec->x0);

        for (resno = 0; resno < tilec->numresolutions; ++resno)
        {
            opj_tcd_resolution_t *res = &tilec->resolutions[resno];

            for (bandno = 0; bandno < res->numbands; ++bandno)
            {
                opj_tcd_band_t* restrict band = &res->bands[bandno];
                OPJ_INT32 bandconst = 8192 * 8192 / ((OPJ_INT32) floor(band->stepsize * 8192));

                OPJ_UINT32 level = tilec->numresolutions - 1 - resno ;//range 0-5
                OPJ_UINT32 orient = band->bandno;

                for (precno = 0; precno < res->pw * res->ph; ++precno)
                {
                    opj_tcd_precinct_t *prc = &band->precincts[precno];

                    for (cblkno = 0; cblkno < prc->cw * prc->ch; ++cblkno)
                    {
                        opj_tcd_cblk_enc_t* cblk = &prc->cblks.enc[cblkno];
                        OPJ_INT32 * restrict datap;
                        OPJ_INT32* restrict tiledp;
                        OPJ_UINT32 cblk_w;
                        OPJ_UINT32 cblk_h;
                        OPJ_UINT32 i, j;

                        OPJ_INT32 x = cblk->x0 - band->x0;
                        OPJ_INT32 y = cblk->y0 - band->y0;
                        if (band->bandno & 1)
                        {
                            opj_tcd_resolution_t *pres = &tilec->resolutions[resno - 1];
                            x += pres->x1 - pres->x0;
                        }
                        if (band->bandno & 2)
                        {
                            opj_tcd_resolution_t *pres = &tilec->resolutions[resno - 1];
                            y += pres->y1 - pres->y0;
                        }

                        if(!opj_t1_allocate_buffers(
                                    t1,
                                    (OPJ_UINT32)(cblk->x1 - cblk->x0),
                                    (OPJ_UINT32)(cblk->y1 - cblk->y0)))
                        {
                            return OPJ_FALSE;
                        }

                        datap=t1->data;
                        cblk_w = t1->w;
                        cblk_h = t1->h;

                        tiledp=&tilec->data[(OPJ_UINT32)y * tile_w + (OPJ_UINT32)x];
                        if (tccp->qmfbid == 1)
                        {
                            for (j = 0; j < cblk_h; ++j)
                            {
                                for (i = 0; i < cblk_w; ++i)
                                {
                                    OPJ_INT32 tmp = tiledp[(j * tile_w) + i];
                                    datap[(j * cblk_w) + i] = tmp << T1_NMSEDEC_FRACBITS;
                                }
                            }
                        }
                        else  		/* if (tccp->qmfbid == 0) */
                        {
#ifdef IS_DEBUG_PRINTF
                            fprintf(stderr,"=====================compno=%d,level=%d,orient=%d,precno=%d,cblkno=%d(%d,%d) ===============================\n",compno, tilec->numresolutions - 1 - resno, band->bandno, precno, cblkno,x,y);
                            fprintf(stderr,"[CodeBlock Size]cblk_h=%d,cblk_w=%d\n",cblk_h,cblk_w);
                            fprintf(stderr,"[Quantization Paramter]step_size=%lf,bandconst=%d\n",band->stepsize,bandconst);
                            fprintf(stderr,"[Layering Paramter]num_vts%d\n",tccp->num_vts);

                            char filename[100];
                            sprintf(filename, "exp_dwt_coeffs_comp%d_level%d_orient%d_codeblock%d", compno,level+1,orient,cblkno,cblk->x0,cblk->y0);
                            fprintf(stderr,"writing DWT coeffs to %s\n",filename);
                            FILE *fp_dwt_coeffs;
                            int min_write_coeff_level=0;
                            if (level>=min_write_coeff_level)
                            {
                                fp_dwt_coeffs=fopen(filename, "wb");
                                fwrite(&cblk_h,sizeof(cblk_h),1,fp_dwt_coeffs);
                                fwrite(&cblk_w,sizeof(cblk_w),1,fp_dwt_coeffs);
                                fwrite(&(cblk->y0),sizeof(cblk->y0),1,fp_dwt_coeffs);
                                fwrite(&(cblk->x0),sizeof(cblk->x0),1,fp_dwt_coeffs);
                            }
#endif
                            for (j = 0; j < cblk_h; ++j)
                            {

                                for (i = 0; i < cblk_w; ++i)
                                {
                                    OPJ_INT32 tmp = tiledp[(j * tile_w) + i];

                                    datap[(j * cblk_w) + i] =
                                        opj_int_fix_mul(
                                            tmp,
                                            bandconst) >> (11 - T1_NMSEDEC_FRACBITS);

#ifdef IS_DEBUG_PRINTF
                                    /*yl: tracking which resno,orient,cblkno,x,y) is the current data mapped to*/

                                    OPJ_UINT32 idx = (OPJ_UINT32)y * tile_w + (OPJ_UINT32)x + (j * tile_w) + i;
                                    if (compno==0)
                                    {
                                        fwrite(&idx,sizeof(OPJ_UINT32),1,fp_dwt_idx); /*idx*/
                                        fwrite(&resno,sizeof(OPJ_UINT32),1,fp_dwt_idx); /*resno*/
                                        fwrite(&bandno,sizeof(OPJ_UINT32),1,fp_dwt_idx); /*bandno*/
                                        fwrite(&cblkno,sizeof(OPJ_UINT32),1,fp_dwt_idx); /*cblkno*/
                                        fwrite(&i,sizeof(OPJ_UINT32),1,fp_dwt_idx); /*x*/
                                        fwrite(&j,sizeof(OPJ_UINT32),1,fp_dwt_idx); /*y*/
                                    }
                                    if (level>=min_write_coeff_level)
                                        fwrite(&tmp, sizeof(tmp),1, fp_dwt_coeffs);
#endif


//                                      fprintf(stderr,"%d,%d\n",tmp,datap[(j * cblk_w) + i]);
                                }
                            }
#ifdef IS_DEBUG_PRINTF

                            if (level>=min_write_coeff_level)
                                fclose(fp_dwt_coeffs);
#endif
                        }
//yl
                        //opj_t1_encode_cblk(
                        if (!tccp->custom_stepsize)
                            opj_t1_encode_cblk(
                                t1,
                                cblk,
                                band->bandno,
                                compno,
                                tilec->numresolutions - 1 - resno,
                                tccp->qmfbid,
                                band->stepsize,
                                tccp->cblksty,
                                tile->numcomps,
                                tile,
                                mct_norms);
                        else
                            yl_opj_t1_encode_cblk(
                                t1,
                                cblk,
                                band->bandno,
                                compno,
                                tilec->numresolutions - 1 - resno,
                                tccp->qmfbid,
                                band->stepsize,
                                tccp->cblksty,
                                tile->numcomps,
                                tile,
                                mct_norms,
                                tccp->lin_mod_slope,
                                tccp->lin_mod_intercept,
                                tccp->num_vts,
                                tiledp);

//yl:end
                    } /* cblkno */
                } /* precno */
            } /* bandno */
        } /* resno  */
    } /* compno  */
#ifdef IS_DEBUG_PRINTF
    fclose(fp_dwt_idx);
#endif
    return OPJ_TRUE;
}



void yl_opj_t1_encode_cblk(opj_t1_t *t1,
                           opj_tcd_cblk_enc_t* cblk,
                           OPJ_UINT32 orient,
                           OPJ_UINT32 compno,
                           OPJ_UINT32 level,
                           OPJ_UINT32 qmfbid,
                           OPJ_FLOAT64 stepsize,
                           OPJ_UINT32 cblksty,
                           OPJ_UINT32 numcomps,
                           opj_tcd_tile_t * tile,
                           const OPJ_FLOAT64 * mct_norms,
                           OPJ_FLOAT64 *lin_mod_slope,
                           OPJ_FLOAT64 *lin_mod_intercept,
                           OPJ_INT32 num_vts,
                           OPJ_INT32* tiledp)
{
#ifdef IS_DEBUG_PRINTF
    fprintf(stderr,"[T1 encoding]:yl_opj_t1_encode_cblk(orient=%d,level=%d)\n",orient,level);
    fprintf(stderr,"allocating memory for decoded sample:size=%d\n",t1->w*t1->h);
#endif // IS_DEBUG_PRINTF

    /*Varible holding decoded data*/
    OPJ_INT32 datap_dec[4096];
    memset((OPJ_INT32*)datap_dec,0,sizeof(OPJ_INT32)*4096 );


    OPJ_FLOAT64 cumwmsedec = 0.0;

    opj_mqc_t *mqc = t1->mqc;	/* MQC component */

    OPJ_UINT32 passno;
    OPJ_INT32 bpno;
    OPJ_UINT32 passtype;
    OPJ_INT32 nmsedec = 0;
    OPJ_UINT32 i,j;
    OPJ_BYTE type = T1_TYPE_MQ;
    OPJ_FLOAT64 tempwmsedec;
    OPJ_INT32 max;

    /*Find the wavelet coefficients with largest magnitude*/

    max=0;
    for (i = 0; i < t1->w * t1->h; ++i)
    {
        OPJ_INT32 tmp = abs(t1->data[i]);
        max = opj_int_max(max, tmp);
    }

    cblk->numbps = max ? (OPJ_UINT32)((opj_int_floorlog2(max) + 1) - T1_NMSEDEC_FRACBITS) : 0;
    if (cblk->numbps < 0) cblk->numbps=0;

    bpno = (OPJ_INT32)(cblk->numbps - 1);
    passtype = 2;

    opj_mqc_resetstates(mqc);
    opj_mqc_setstate(mqc, T1_CTXNO_UNI, 0, 46);
    opj_mqc_setstate(mqc, T1_CTXNO_AGG, 0, 3);
    opj_mqc_setstate(mqc, T1_CTXNO_ZC, 0, 4);
    opj_mqc_init_enc(mqc, cblk->data);



#ifdef IS_DEBUG_PRINTF
    fprintf(stderr,"max of code block is %d,cblk->numbps=%d\n",max,cblk->numbps);
#endif


    /*Calculate the following stats of the codeblock(used later to decide jnd_threshold and masking coefficients*/
    /*variance(l2-norm)*/
    OPJ_FLOAT64 sample_var=0;
    /*mean*/
    OPJ_FLOAT64 sample_mean2=0;
    /*mean of magnitudes( l1-norm)*/
    OPJ_FLOAT64 sample_mean=0;
    /*maximum (infinity norm)*/
    OPJ_FLOAT64 max_coeff=-1000.0;

    OPJ_FLOAT64 val_quant;
    OPJ_UINT32 num_coeffs=t1->w * t1->h;
    for (i = 0; i < t1->w; ++i)
    {
        for (j = 0; j < t1->h; ++j)
        {
            //  val_quant=(OPJ_FLOAT64) (opj_int_abs(t1->data[(j * t1->w) + i])>>T1_NMSEDEC_FRACBITS);
            val_quant=(OPJ_FLOAT64) (opj_int_abs(t1->data[(j * t1->w) + i])) / pow(2.0,(OPJ_FLOAT64)T1_NMSEDEC_FRACBITS);
            if (val_quant > max_coeff)
                max_coeff = val_quant;
            sample_mean2+=val_quant*stepsize;
            if (t1->data[(j * t1->w) + i] <0)
                val_quant=-val_quant;
            val_quant*=stepsize;// /64.0*256.0; //Range: (-128,128)
            sample_mean+=val_quant;
            sample_var+=val_quant*val_quant;
        }
    }
    sample_mean=sample_mean/num_coeffs;
    sample_var=sample_var/num_coeffs-sample_mean*sample_mean;
    sample_mean2=sample_mean2/num_coeffs;


    /*calculating threshold from linear model*/
    OPJ_FLOAT64 jnd_threshold[100];
    /*Threshold for determine wether or not distortion should be multiplied by 2*/
    OPJ_FLOAT64 qthresh[100];
    OPJ_UINT32 i_vt=0; /*VT index,for single layer alwasy 0*/
    OPJ_UINT32 i_band;
    for (; i_vt < num_vts; i_vt++)
    {
        //For Luminance component, visibility threshold is decided by the linear model
        /*
        if (compno==0 && level<5)
            jnd_threshold[i_vt] = lin_mod_slope[level*3 + orient-1+i_vt*16]*sample_var+lin_mod_intercept[level*3 + orient-1+i_vt*48];
        //For Chrominance, it's constant that is  independent of the sample variance
        else if (level==5)  //Luminance/ Chrominance LL5
            jnd_threshold[i_vt] = lin_mod_intercept[level*3 + orient+i_vt*48];
        else     //Chrominance LH/HL/HH
            jnd_threshold[i_vt] = lin_mod_intercept[level*3 + orient-1+i_vt*48];
        */
        i_band = (orient == 0) ? 0: ( ( 4-level) * 3 +orient );
        if (compno==0 && level<5)
            jnd_threshold[i_vt] = lin_mod_slope[i_band+i_vt*16]*sample_var+lin_mod_intercept[i_band+i_vt*48];
        else  //Luminance/ Chrominance LL5 or chrominance
            jnd_threshold[i_vt] = lin_mod_intercept[i_band+i_vt*48];



        qthresh[i_vt]=jnd_threshold[i_vt]/stepsize;
    }

#ifdef IS_DEBUG_PRINTF
    if (compno==0)
    {
        fprintf(stderr,"(compno=%d,level=%d,orient=%d),sample_var=%lf\n,sample_mean2=%lf,max_coeff=%lf",compno,level,orient,sample_var,sample_mean2,max_coeff);
        for (i_vt=0 ; i_vt<num_vts; i_vt++)
            fprintf(stderr,",jnd_threshold_no_masking[%d]=%lf=%lf sigma^2 +%lf",i_vt,jnd_threshold[i_vt],lin_mod_slope[level*3 + orient-1+i_vt*16],lin_mod_intercept[level*3 + orient-1+i_vt*48]);
        fprintf(stderr,"\n");
    }
    else
    {
        fprintf(stderr,"(compno=%d,level=%d,orient=%d),jnd_threshold_no_mask=%lf,sample_var=%lf\n,sample_mean2=%lf,max_coeff=%lf",compno,level,orient,jnd_threshold,sample_var,sample_mean2,max_coeff);
        for (i_vt=0 ; i_vt<num_vts; i_vt++)
            fprintf(stderr,",jnd_threshold_no_masking[%d]=%lf=%lf sigma^2 +%lf",i_vt,jnd_threshold[i_vt],0.0,lin_mod_intercept[level*3 + orient-1+i_vt*48]);
        fprintf(stderr,"\n");
    }
    //yl:end
#endif

    //yl:masking
#if IS_USING_MASKING
    OPJ_FLOAT64 sm[4096];
    OPJ_FLOAT64 nm[4096];
    memset((OPJ_FLOAT64*)sm,0,sizeof(OPJ_FLOAT64)*4096 );
    memset((OPJ_FLOAT64*)nm,0,sizeof(OPJ_FLOAT64)*4096 );
    //yl:end
    OPJ_FLOAT64 sm_weight = 0.8;
    OPJ_FLOAT64 sm_exp = 0.35;
    for (i=0; i<num_coeffs; i++)
    {
        //val_quant=(OPJ_FLOAT64) (opj_int_abs(t1->data[(j * t1->w) + i])>>T1_NMSEDEC_FRACBITS);
        val_quant=(OPJ_FLOAT64) (opj_int_abs(t1->data[i])) / pow(2.0,(OPJ_FLOAT64)T1_NMSEDEC_FRACBITS); //fix the variance descrapancy
        val_quant*=stepsize;
        sm[i] = sm_weight*pow(val_quant/(sample_mean2+0.00001),sm_exp); //equ(8)
        if (sm[i]<1.0)
            sm[i]=1.0;
    }

    //yl:end

    //yl:texture-masking
    OPJ_INT32 winsize = 32;
    OPJ_INT32 winsizeb = (int)((OPJ_FLOAT64)winsize/pow(2.0,level+1));
    OPJ_FLOAT64 nm_weight = 0.35;
    OPJ_FLOAT64 nm_exp = 0.4;
    OPJ_FLOAT64 masking_beta = 0.8;
    OPJ_FLOAT64 masking_factor=1.0;
    fprintf(stderr,"coding passes.....\n");
//yl:end

#endif
    OPJ_INT32 last_pass_vt_bin=-2;


    /*newest added varialbes for slope logic*/
    OPJ_FLOAT64 minmaxdist=1000.0;
    OPJ_FLOAT64 firstlow[100];
    OPJ_UINT32 minmaxdistidx;
    OPJ_UINT32 firstlowidx[100];
    OPJ_BOOL below_threshold[100];
    for (i_vt=0; i_vt < num_vts; i_vt++)
    {
        below_threshold[i_vt] = (0==1);
        firstlowidx[i_vt] = 1000;
    }




    for (passno = 0; bpno >= 0; ++passno)
    {
        opj_tcd_pass_t *pass = &cblk->passes[passno];
        OPJ_UINT32 correction = 3;
        type = ((bpno < ((OPJ_INT32) (cblk->numbps) - 4)) && (passtype < 2) && (cblksty & J2K_CCP_CBLKSTY_LAZY)) ? T1_TYPE_RAW : T1_TYPE_MQ;

        switch (passtype)
        {
        case 0:
            yl_opj_t1_enc_sigpass(t1,datap_dec, bpno, orient, &nmsedec, type, cblksty);
            break;
        case 1:
            yl_opj_t1_enc_refpass(t1,datap_dec, bpno, &nmsedec, type, cblksty);
            break;
        case 2:
            yl_opj_t1_enc_clnpass(t1,datap_dec, bpno, orient, &nmsedec, cblksty);
            /* code switch SEGMARK (i.e. SEGSYM) */
            if (cblksty & J2K_CCP_CBLKSTY_SEGSYM)
                opj_mqc_segmark_enc(mqc);
            break;
        }



        opj_tcd_tilecomp_t* tilec = &tile->comps[compno];
        OPJ_UINT32 tile_w = (OPJ_UINT32)(tilec->x1 - tilec->x0);

        OPJ_UINT32 cblk_w = t1->w;
        OPJ_UINT32 cblk_h = t1->h;

        /*find the largest distortion in the codeblock*/
        OPJ_FLOAT64 max_dist[100];
        OPJ_FLOAT64 max_diff[100];
        for (i_vt=0; i_vt<100; i_vt++)
            max_diff[i_vt]=-1.0;
        OPJ_FLOAT64 dist;
        OPJ_FLOAT64 diff;
        OPJ_FLOAT64 max_diff_dist_gain[100];
        OPJ_FLOAT64 x[100];
        OPJ_FLOAT64 x_dec[100];
        //OPJ_INT32 y,y_dec;
        OPJ_FLOAT64 y,y_dec; //Feng's newest diortion calculation uses floating point
        /*distortion gain (or not) depending wether coeffciients is larger than q_thres*/
        OPJ_FLOAT64 dist_gain[100]; //2.0 or 1.0


        /*find the largest distortion in the codeblock*/
        for (i = 0; i < t1->w * t1->h; i++ )
        {
            //diff=opj_int_abs(t1->data [i] - datap_dec[ i]);
            //y = opj_int_abs(t1->data [i]) >>T1_NMSEDEC_FRACBITS ;
            //y_dec = opj_int_abs(datap_dec[i]) >>T1_NMSEDEC_FRACBITS ;

            y = (OPJ_FLOAT64) (opj_int_abs(t1->data[i])) / pow(2.0,(OPJ_FLOAT64)T1_NMSEDEC_FRACBITS); /* coeff. without bit plane trucation*/
            y_dec = (OPJ_FLOAT64) (opj_int_abs(datap_dec[i])) / pow(2.0,(OPJ_FLOAT64)T1_NMSEDEC_FRACBITS);/*coeff. truncated at current coding pass*/

            for (i_vt =0; i_vt<num_vts; i_vt++)
            {
                if (y >=qthresh[i_vt] && y_dec>=qthresh[i_vt])
                {
                    dist_gain[i_vt]=2.0;
                }
                else
                {
                    dist_gain[i_vt]=1.0;
                }
            }
            if (t1->data[i]<0)
                y=-y;
            if (datap_dec[i]<0)
                y_dec=-y_dec;
            diff = y>y_dec?y-y_dec:y_dec-y;


            for (i_vt =0; i_vt<num_vts; i_vt++)
            {
                if  (dist_gain[i_vt]*diff > max_diff[i_vt])
                {
                    max_diff[i_vt]=dist_gain[i_vt]*diff;
                    x_dec[i_vt]=y_dec;
                    x[i_vt]=y;
                    max_diff_dist_gain[i_vt] = dist_gain[i_vt];
                }
            }

        }
        // }

        //max_dist=(OPJ_FLOAT64) (opj_int_abs(max_diff)>>T1_NMSEDEC_FRACBITS)*stepsize;//Numerical Instability!!! when it is 1x stepsize
        // max_dist=(OPJ_FLOAT64) (opj_int_abs(max_diff))*stepsize;//Numerical Instability!!! when it is 1x stepsize

        for (i_vt =0; i_vt<num_vts; i_vt++)
            max_dist[i_vt]=max_diff[i_vt]*stepsize;
        //yl: Solve numerical instability introduced by using floating point operation
        //if (level==5 || compno!=0 )
        //if ((opj_int_abs(max_diff)) ==1){
        //        if ((opj_int_abs(max_diff)>>T1_NMSEDEC_FRACBITS) ==1){
        //max_dist = jnd_threshold[num_vts-1]; //assumes jnd_threshold[num_vts-1]==stepsize mathematically
        //fprintf(stderr,"rounding to integer jnd\n");
        //}

        /*see if the threshold needs to be adjusted because of : (1) masking effect; (2) multi-layer codestream*/
        OPJ_FLOAT64 adjusted_threshold[100];
        for (i_vt=0; i_vt<num_vts; i_vt++)
        {
            adjusted_threshold[i_vt] = jnd_threshold[i_vt];
            if (i_vt > 0) //ONLY FOR MULTI-LAYER:Maker sure that threshold is in desceding order. if not trust the
            {
                if ( adjusted_threshold[i_vt] > adjusted_threshold[i_vt-1])
                {
                    adjusted_threshold[i_vt] = adjusted_threshold[i_vt-1]; //Between VT_i,and VT_i+1, trust the smaller one.
                    fprintf(stderr,"[Warning]:layer%d's threshold is higher than layer%d's\n",i_vt-1,i_vt);
                }
            }
        }
#if IS_USING_MASKING
        //calculating texture-masking
        OPJ_INT32 m,n;
        OPJ_FLOAT64 sample_val,nm_val;

        if (level <= 2  && (compno==0))
        {
            for (i = 0; i <  t1->w; i += winsizeb)
            {
                for ( j = 0; j < t1->h; j += winsizeb)
                {
                    sample_var = 0;
                    sample_mean = 0;
                    //calculate local variance
                    for (m = i; m < i+winsizeb; m++)  //m,y,i - (row)
                    {
                        for (n = j; n < j+winsizeb; n++) //n,x,j - (col)
                        {
                            //sample_val = (OPJ_FLOAT64) (opj_int_abs(datap_dec[n * t1->w + m])>>T1_NMSEDEC_FRACBITS);
                            sample_val=(OPJ_FLOAT64) (opj_int_abs(datap_dec[n*t1->w+m])) / pow(2.0,(OPJ_FLOAT64)T1_NMSEDEC_FRACBITS);
                            sample_val=sample_val*stepsize;
                            //  (OPJ_FLOAT64)(dec_samples[m*block->size.get_x()+n]&0x7fffffff)/(OPJ_FLOAT64)(1<<(31-block->k_max));
                            if (datap_dec[n * t1->w + m] < 0)
                                sample_val = -sample_val;
                            sample_var += sample_val*sample_val;
                            sample_mean += sample_val;
                        }
                    }
                    sample_mean = 1.0/(OPJ_FLOAT64)(winsizeb*winsizeb)*sample_mean;
                    sample_var = 1.0/(OPJ_FLOAT64)(winsizeb*winsizeb)*sample_var - sample_mean*sample_mean;

                    //sample_var *= stepsize;
                    nm_val = nm_weight*pow(sample_var,nm_exp);
                    if (nm_val < 1.0)
                        nm_val = 1.0;

                    //allocate texture masking
                    for ( m = i; m < i+winsizeb; m++)
                    {
                        for ( n = j; n < j+winsizeb; n++)
                        {
                            nm[n * t1->w+ m] = nm_val;
                        }
                    }

                }
            }

            sample_mean = 0;
            OPJ_FLOAT64 masking;
            //minkowski mean
            for ( i = 0; i < num_coeffs; i++)
            {
                //if(level==2 && orient==1)
                // fprintf(stderr,"sm[%d]=%lf,nm[%d]=%lf,yhat=%f,y=%f\n",i,sm[i],i,nm[i],(OPJ_FLOAT64) (opj_int_abs(datap_dec[i])) / pow(2.0,(OPJ_FLOAT64)T1_NMSEDEC_FRACBITS),(OPJ_FLOAT64)(opj_int_abs(t1->data[i])) / pow(2.0,(OPJ_FLOAT64)T1_NMSEDEC_FRACBITS));

                masking = sm[i]*nm[i];
                //	  fprintf(stderr,"sm=%lf,nm=%lf\n",sm[i],nm[i]);
                sample_mean += pow(masking,masking_beta);
            }
            masking_factor = pow(1.0/(OPJ_FLOAT64)num_coeffs*sample_mean,1.0/masking_beta);

        }

        //yl:end

        for (i_vt=0; i_vt<num_vts; i_vt++)
        {
            if (level==5 && orient==0)
                adjusted_threshold[i_vt] = masking_factor*jnd_threshold[i_vt];
            else
            {
                OPJ_FLOAT64 b_tmp = lin_mod_intercept[level*3+orient-1+ i_vt * 48];
                adjusted_threshold[i_vt] = ( masking_factor*jnd_threshold[i_vt]< ( 6.0 * b_tmp) ? masking_factor*jnd_threshold[i_vt] : (6.0*b_tmp) );
            }
        }

#endif // IS_USING_MASKING

        /*newest switch for slope logic*/
        /*catch the coding pass where the actual minimal distortin occurred(not necessarily the last coding pass)*/
        if((max_dist[num_vts]<minmaxdist)||(passno==0))
        {
            minmaxdist=max_dist[num_vts];
            minmaxdistidx=passno;
#ifdef IS_DEBUG_PRINTF
            fprintf(stderr,"===============min { max_dist updated }==================\n");
#endif
        }
        /*catch the first coding pass where the max distortion is below the threshold*/
        for (i_vt=0; i_vt < num_vts; i_vt++)
        {
            if((!below_threshold[i_vt] )&&(max_dist[i_vt]<=adjusted_threshold[i_vt]))
            {
                firstlow[i_vt]=max_dist[i_vt];
                firstlowidx[i_vt]=passno;
                below_threshold[i_vt]=(1==1);
            }
        }


#ifdef IS_DEBUG_PRINTF
        for (i_vt=0; i_vt<num_vts; i_vt++)
        {
            fprintf(stderr,"[pass%d],max_dist=%lf=|%f-%f|*%f*%f,masking_factor=%lf",passno,max_dist[i_vt],x[i_vt],x_dec[i_vt],max_diff_dist_gain[i_vt],stepsize,masking_factor);
            // fprintf(stderr,",coeff[344]=%lf,dec[344]=%lf",(OPJ_FLOAT64) (opj_int_abs(t1->data[344])) / pow(2.0,(OPJ_FLOAT64)T1_NMSEDEC_FRACBITS),(OPJ_FLOAT64) (opj_int_abs(datap_dec[344])) / pow(2.0,(OPJ_FLOAT64)T1_NMSEDEC_FRACBITS));
            // fprintf(stderr,",coeff[345]=%lf,dec[345]=%lf",(OPJ_FLOAT64) (opj_int_abs(t1->data[343])) / pow(2.0,(OPJ_FLOAT64)T1_NMSEDEC_FRACBITS),(OPJ_FLOAT64) (opj_int_abs(datap_dec[343])) / pow(2.0,(OPJ_FLOAT64)T1_NMSEDEC_FRACBITS));
            fprintf(stderr,",adjusted_threshold[%d]=%lf\n",i_vt,adjusted_threshold[i_vt]);
        }
#endif
        //if (max_dist<=adjusted_threshold[num_vts-1]) //Distortion <= VT1(stepsize)
        if (below_threshold[num_vts-1] && firstlowidx[num_vts-1] !=passno)
        {
            /*below threshold, but not the first time. The coding pass for first time crossing the threshld is always coded to be coservative*/
#ifdef IS_DEBUG_PRINTF
            fprintf(stderr,"DISCARD BIT PLANE");
#endif
            pass->vt_bin=-1; /*this pass goes to layer -1, i.e. discarding this coding pass*/
        }
        else
        {
            OPJ_INT32 i_tmp=0;
            //if (max_dist>adjusted_threshold[0]) //Distortion larger than the lagerste VT
            if (!below_threshold[0]|| (below_threshold[0] && passno==firstlowidx[0]) )
                pass->vt_bin=i_tmp; /*encode the coding pass if :(1) max dist. not smaller than threshold or (2) first time max dist smaller than threshold
            else/*only for multi-layer*/
            {
                for (; i_tmp<num_vts-1; i_tmp++)
                {
                    //    if ( max_dist > adjusted_threshold[i_tmp+1] && max_dist <= adjusted_threshold[i_tmp] ){
                    if ( ( (!below_threshold[i_tmp+1]) || (below_threshold[i_tmp+1] && passno == firstlowidx[i_tmp+1]) )
                            && (below_threshold[i_tmp]   && passno != firstlowidx[i_tmp]    ) )
                    {
                        pass->vt_bin = i_tmp+1;
                        break;
                    }
                }
                if (i_tmp == num_vts-2)
                {
                    if (below_threshold[num_vts-1] && firstlowidx[num_vts-1]==passno)
                        pass->vt_bin = num_vts-1; //newest fix, was num_vts
                }
            }
        }
#ifdef IS_DEBUG_PRINTF
        fprintf(stderr,"=====>layno=%d\n",pass->vt_bin);
#endif


        /* fixed_quality */
        tempwmsedec = opj_t1_getwmsedec(nmsedec, compno, level, orient, bpno, qmfbid, stepsize, numcomps,mct_norms) ;
        cumwmsedec += tempwmsedec;
        tile->distotile += tempwmsedec;

        /* Code switch "RESTART" (i.e. TERMALL) */
        if ((cblksty & J2K_CCP_CBLKSTY_TERMALL)	&& !((passtype == 2) && (bpno - 1 < 0)))
        {
            if (type == T1_TYPE_RAW)
            {
                opj_mqc_flush(mqc);
                correction = 1;
                /* correction = mqc_bypass_flush_enc(); */
            }
            else  			/* correction = mqc_restart_enc(); */
            {
                opj_mqc_flush(mqc);
                correction = 1;
            }
            pass->term = 1;
        }
        else
        {
            if (((bpno < ((OPJ_INT32) (cblk->numbps) - 4) && (passtype > 0))
                    || ((bpno == ((OPJ_INT32)cblk->numbps - 4)) && (passtype == 2))) && (cblksty & J2K_CCP_CBLKSTY_LAZY))
            {
                if (type == T1_TYPE_RAW)
                {
                    opj_mqc_flush(mqc);
                    correction = 1;
                    /* correction = mqc_bypass_flush_enc(); */
                }
                else  		/* correction = mqc_restart_enc(); */
                {
                    opj_mqc_flush(mqc);
                    correction = 1;
                }
                pass->term = 1;
            }
            else
            {
                pass->term = 0;
            }
        }

        if (++passtype == 3)
        {
            passtype = 0;
            bpno--;
        }

        if (pass->term && bpno > 0)
        {
            type = ((bpno < ((OPJ_INT32) (cblk->numbps) - 4)) && (passtype < 2) && (cblksty & J2K_CCP_CBLKSTY_LAZY)) ? T1_TYPE_RAW : T1_TYPE_MQ;
            if (type == T1_TYPE_RAW)
                opj_mqc_bypass_init_enc(mqc);
            else
                opj_mqc_restart_init_enc(mqc);
        }

        pass->distortiondec = cumwmsedec;
        pass->rate = opj_mqc_numbytes(mqc) + correction;	/* FIXME */

        /* Code-switch "RESET" */
        if (cblksty & J2K_CCP_CBLKSTY_RESET)
            opj_mqc_reset_enc(mqc);
    } /*end of coding pass loop*/



    /* Code switch "ERTERM" (i.e. PTERM) */
    if (cblksty & J2K_CCP_CBLKSTY_PTERM)
        opj_mqc_erterm_enc(mqc);
    else /* Default coding */ if (!(cblksty & J2K_CCP_CBLKSTY_LAZY))
        opj_mqc_flush(mqc);

    cblk->totalpasses = passno;

    for (passno = 0; passno<cblk->totalpasses; passno++)
    {
        opj_tcd_pass_t *pass = &cblk->passes[passno];
        if (pass->rate > opj_mqc_numbytes(mqc))
            pass->rate = opj_mqc_numbytes(mqc);
        /*Preventing generation of FF as last data byte of a pass*/
        if((pass->rate>1) && (cblk->data[pass->rate - 1] == 0xFF))
        {
            pass->rate--;
        }
        pass->len = pass->rate - (passno == 0 ? 0 : cblk->passes[passno - 1].rate);
    }

    /*newest slope logic */
    /* there might be repeatation within this part that is already done before*/
    for (passno = 0; passno<cblk->totalpasses; passno++)
    {
        opj_tcd_pass_t *pass = &cblk->passes[passno];
        if (max_coeff < jnd_threshold[num_vts-1] && pass->vt_bin == num_vts)
            pass->vt_bin = -1;     /*if the max coefficients is smaller than the threshold, discard ALL the coding passes*/
        for (i_vt = 0; i_vt < num_vts; i_vt ++)
        {
            if (!below_threshold[i_vt] && passno<=minmaxdistidx)
            {
                pass->vt_bin=i_vt; /* ONly encode pass that is before the pass where actural minimal distrotion occurs*/
                break; /*only for multi-layer*/
            }
        }
#ifdef IS_DEBUG_PRINTF
        fprintf(stderr,"pass[%d](rate=%d,mse=%f)=====>layno=%d\n",passno,pass->rate,pass->distortiondec,pass->vt_bin);
#endif
    }

}






/** mod fixed_quality */
void opj_t1_encode_cblk(opj_t1_t *t1,
                        opj_tcd_cblk_enc_t* cblk,
                        OPJ_UINT32 orient,
                        OPJ_UINT32 compno,
                        OPJ_UINT32 level,
                        OPJ_UINT32 qmfbid,
                        OPJ_FLOAT64 stepsize,
                        OPJ_UINT32 cblksty,
                        OPJ_UINT32 numcomps,
                        opj_tcd_tile_t * tile,
                        const OPJ_FLOAT64 * mct_norms)
{
    fprintf(stderr,"default tier one coding....\n");
    OPJ_FLOAT64 cumwmsedec = 0.0;

    opj_mqc_t *mqc = t1->mqc;	/* MQC component */

    OPJ_UINT32 passno;
    OPJ_INT32 bpno;
    OPJ_UINT32 passtype;
    OPJ_INT32 nmsedec = 0;
    OPJ_INT32 max;
    OPJ_UINT32 i;
    OPJ_BYTE type = T1_TYPE_MQ;
    OPJ_FLOAT64 tempwmsedec;

    max = 0;
    for (i = 0; i < t1->w * t1->h; ++i)
    {
        OPJ_INT32 tmp = abs(t1->data[i]);
        max = opj_int_max(max, tmp);
    }

    cblk->numbps = max ? (OPJ_UINT32)((opj_int_floorlog2(max) + 1) - T1_NMSEDEC_FRACBITS) : 0;
    fprintf(stderr,"cblk->numbps=%d\n",cblk->numbps);

    bpno = (OPJ_INT32)(cblk->numbps - 1);
    passtype = 2;

    opj_mqc_resetstates(mqc);
    opj_mqc_setstate(mqc, T1_CTXNO_UNI, 0, 46);
    opj_mqc_setstate(mqc, T1_CTXNO_AGG, 0, 3);
    opj_mqc_setstate(mqc, T1_CTXNO_ZC, 0, 4);
    opj_mqc_init_enc(mqc, cblk->data);




    for (passno = 0; bpno >= 0; ++passno)
    {
        opj_tcd_pass_t *pass = &cblk->passes[passno];
        OPJ_UINT32 correction = 3;
        type = ((bpno < ((OPJ_INT32) (cblk->numbps) - 4)) && (passtype < 2) && (cblksty & J2K_CCP_CBLKSTY_LAZY)) ? T1_TYPE_RAW : T1_TYPE_MQ;

        switch (passtype)
        {
        case 0:
            opj_t1_enc_sigpass(t1, bpno, orient, &nmsedec, type, cblksty);
            break;
        case 1:
            opj_t1_enc_refpass(t1, bpno, &nmsedec, type, cblksty);
            break;
        case 2:
            opj_t1_enc_clnpass(t1, bpno, orient, &nmsedec, cblksty);
            /* code switch SEGMARK (i.e. SEGSYM) */
            if (cblksty & J2K_CCP_CBLKSTY_SEGSYM)
                opj_mqc_segmark_enc(mqc);
            break;
        }




        /* fixed_quality */
        tempwmsedec = opj_t1_getwmsedec(nmsedec, compno, level, orient, bpno, qmfbid, stepsize, numcomps,mct_norms) ;
        cumwmsedec += tempwmsedec;
        tile->distotile += tempwmsedec;

        /* Code switch "RESTART" (i.e. TERMALL) */
        if ((cblksty & J2K_CCP_CBLKSTY_TERMALL)	&& !((passtype == 2) && (bpno - 1 < 0)))
        {
            if (type == T1_TYPE_RAW)
            {
                opj_mqc_flush(mqc);
                correction = 1;
                /* correction = mqc_bypass_flush_enc(); */
            }
            else  			/* correction = mqc_restart_enc(); */
            {
                opj_mqc_flush(mqc);
                correction = 1;
            }
            pass->term = 1;
        }
        else
        {
            if (((bpno < ((OPJ_INT32) (cblk->numbps) - 4) && (passtype > 0))
                    || ((bpno == ((OPJ_INT32)cblk->numbps - 4)) && (passtype == 2))) && (cblksty & J2K_CCP_CBLKSTY_LAZY))
            {
                if (type == T1_TYPE_RAW)
                {
                    opj_mqc_flush(mqc);
                    correction = 1;
                    /* correction = mqc_bypass_flush_enc(); */
                }
                else  		/* correction = mqc_restart_enc(); */
                {
                    opj_mqc_flush(mqc);
                    correction = 1;
                }
                pass->term = 1;
            }
            else
            {
                pass->term = 0;
            }
        }

        if (++passtype == 3)
        {
            passtype = 0;
            bpno--;
        }

        if (pass->term && bpno > 0)
        {
            type = ((bpno < ((OPJ_INT32) (cblk->numbps) - 4)) && (passtype < 2) && (cblksty & J2K_CCP_CBLKSTY_LAZY)) ? T1_TYPE_RAW : T1_TYPE_MQ;
            if (type == T1_TYPE_RAW)
                opj_mqc_bypass_init_enc(mqc);
            else
                opj_mqc_restart_init_enc(mqc);
        }

        pass->distortiondec = cumwmsedec;
        pass->rate = opj_mqc_numbytes(mqc) + correction;	/* FIXME */

        /* Code-switch "RESET" */
        if (cblksty & J2K_CCP_CBLKSTY_RESET)
            opj_mqc_reset_enc(mqc);
    }

    /* Code switch "ERTERM" (i.e. PTERM) */
    if (cblksty & J2K_CCP_CBLKSTY_PTERM)
        opj_mqc_erterm_enc(mqc);
    else /* Default coding */ if (!(cblksty & J2K_CCP_CBLKSTY_LAZY))
        opj_mqc_flush(mqc);

    cblk->totalpasses = passno;

    for (passno = 0; passno<cblk->totalpasses; passno++)
    {
        opj_tcd_pass_t *pass = &cblk->passes[passno];
        if (pass->rate > opj_mqc_numbytes(mqc))
            pass->rate = opj_mqc_numbytes(mqc);
        /*Preventing generation of FF as last data byte of a pass*/
        if((pass->rate>1) && (cblk->data[pass->rate - 1] == 0xFF))
        {
            pass->rate--;
        }
        pass->len = pass->rate - (passno == 0 ? 0 : cblk->passes[passno - 1].rate);
    }
}

#if 0
void opj_t1_dec_refpass_step(   opj_t1_t *t1,
                                opj_flag_t *flagsp,
                                OPJ_INT32 *datap,
                                OPJ_INT32 poshalf,
                                OPJ_INT32 neghalf,
                                OPJ_BYTE type,
                                OPJ_UINT32 vsc)
{
    OPJ_INT32  t;
    OPJ_UINT32 v,flag;

    opj_mqc_t *mqc = t1->mqc;	/* MQC component */
    opj_raw_t *raw = t1->raw;	/* RAW component */

    flag = vsc ? ((*flagsp) & (~(T1_SIG_S | T1_SIG_SE | T1_SIG_SW | T1_SGN_S))) : (*flagsp);
    if ((flag & (T1_SIG | T1_VISIT)) == T1_SIG)
    {
        opj_mqc_setcurctx(mqc, opj_t1_getctxno_mag(flag));	/* ESSAI */
        if (type == T1_TYPE_RAW)
        {
            v = opj_raw_decode(raw);
        }
        else
        {
            v = opj_mqc_decode(mqc);
        }
        t = v ? poshalf : neghalf;
        *datap += *datap < 0 ? -t : t;
        *flagsp |= T1_REFINE;
    }
}				/* VSC and  BYPASS by Antonin  */
#endif



#if 0
void opj_t1_dec_sigpass_step(   opj_t1_t *t1,
                                opj_flag_t *flagsp,
                                OPJ_INT32 *datap,
                                OPJ_UINT32 orient,
                                OPJ_INT32 oneplushalf,
                                OPJ_BYTE type,
                                OPJ_UINT32 vsc)
{
    OPJ_UINT32 v, flag;

    opj_raw_t *raw = t1->raw;	/* RAW component */
    opj_mqc_t *mqc = t1->mqc;	/* MQC component */

    flag = vsc ? ((*flagsp) & (~(T1_SIG_S | T1_SIG_SE | T1_SIG_SW | T1_SGN_S))) : (*flagsp);
    if ((flag & T1_SIG_OTH) && !(flag & (T1_SIG | T1_VISIT)))
    {
        if (type == T1_TYPE_RAW)
        {
            if (opj_raw_decode(raw))
            {
                v = opj_raw_decode(raw);	/* ESSAI */
                *datap = v ? -oneplushalf : oneplushalf;
                opj_t1_updateflags(flagsp, v, t1->flags_stride);
            }
        }
        else
        {
            opj_mqc_setcurctx(mqc, opj_t1_getctxno_zc(flag, orient));
            if (opj_mqc_decode(mqc))
            {
                opj_mqc_setcurctx(mqc, opj_t1_getctxno_sc(flag));
                v = opj_mqc_decode(mqc) ^ opj_t1_getspb(flag);
                *datap = v ? -oneplushalf : oneplushalf;
                opj_t1_updateflags(flagsp, v, t1->flags_stride);
            }
        }
        *flagsp |= T1_VISIT;
    }
}				/* VSC and  BYPASS by Antonin */
#endif
