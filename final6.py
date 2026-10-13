# -*- coding: utf-8 -*-
"""Interface Pydroid/Tkinter para o núcleo ZAB v2 (quasipartículas, carga/descarga)."""
import math, os, time, threading, webbrowser
import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox
import tkinter.font as tkfont
# ===== NÚCLEO MATEMÁTICO INCORPORADO (arquivo único) =====
# -*- coding: utf-8 -*-
"""
ZAB neutra ZnCl2-NH4Cl — núcleo v2 (numpy puro).
Base: Clark, Latz e Horstmann (2017) + SI. Aproximações explícitas no fim do arquivo (PREMISSAS).

Estado por volume (mol/m3 de eletrólito): quasipartículas Q = [Zn~, NH3~, NH4~]  (Eqs. 21-23 do artigo)
    N = NH3~ + NH4~ ;  Cl = 2 Zn~ + NH4~   (eletroneutralidade, Eqs. 24-26)
Fontes (Tabela 5): Zn -> Zn~ + 2e- (eletrodo de Zn) ; 2 NH4~ + 2e- -> 2 NH3~ (GDE, ORR; inverso na OER).
Transporte: Nernst-Planck por espécie, potencial resolvido por corrente imposta, difusão implícita (bloco 3x3).
Precipitação: equilíbrio local (SI=0) com 4 sólidos; dissolve na carga.
"""

F = 96485.33212
R = 8.314462618
T = 298.15
RTF = R*T/F
LN10 = math.log(10.0)

# ---------------- espécies (mesma tabela/ordem do final5, 16 aquosas) ----------------
# (nome, M g/mol, z, D m2/s, lambda S m2/equiv, Vbar m3/mol)
SPECIES = [
    ("ZnCl+", 100.8,  1, 1.10e-9, 35.0e-4, -1.28e-6),
    ("ZnCl2",  136.3, 0, 1.10e-9,  0.0,    24.82e-6),
    ("ZnCl3-",171.7, -1, 1.10e-9, 35.0e-4,  53.9e-6),
    ("ZnCl4--",207.2,-2, 1.10e-9, 55.0e-4,  81.0e-6),
    ("Zn(NH3)2+",82.4, 2, 1.12e-9, 84.2e-4,  7.5e-6),
    ("Zn(NH3)2^2+",99.4,2,1.14e-9, 85.7e-4, 39.2e-6),
    ("Zn(NH3)3^2+",116.5,2,1.09e-9,81.5e-4, 70.9e-6),
    ("Zn(NH3)4^2+",133.5,2,0.83e-9,62.7e-4,102.6e-6),
    ("ZnCl3(NH3)-",188.8,-1,1.10e-9,41.3e-4, 86.4e-6),
    ("ZnCl(NH3)3+",151.9,1,1.10e-9,41.3e-4,  97.2e-6),
    ("H+",1.0,1,9.31e-9,349.8e-4,0.0),
    ("OH-",17.0,-1,5.26e-9,197.6e-4,-4.18e-6),
    ("NH4+",18.0,1,1.95e-9,73.4e-4,18.13e-6),
    ("NH3",17.0,0,1.94e-9,0.0,24.43e-6),
    ("Cl-",35.5,-1,2.03e-9,76.34e-4,17.79e-6),
    ("Zn2+",65.4,2,0.71e-9,53.0e-4,-24.3e-6),
]
NS = len(SPECIES)
IDX = {s[0]: i for i, s in enumerate(SPECIES)}
IZN, ICL, INH3, INH4, IH, IOH = IDX["Zn2+"], IDX["Cl-"], IDX["NH3"], IDX["NH4+"], IDX["H+"], IDX["OH-"]
Z_CH = np.array([s[2] for s in SPECIES], float)
D_SP = np.array([s[3] for s in SPECIES], float)
ZD = Z_CH*D_SP
Z2D = Z_CH**2*D_SP

# 12 log10(beta) da Tabela 1
LOG_BETA = [0.10, 0.06, 0.10, 0.30, 2.38, 4.88, 7.43, 9.65, 3.15, 7.90, 9.80, 13.96]
# expoentes de [Zn2+], [Cl-], [NH3], [H+] e log10(beta) por espécie
_EXP = {
    "ZnCl+": ((1,1,0,0), LOG_BETA[0]), "ZnCl2": ((1,2,0,0), LOG_BETA[1]),
    "ZnCl3-": ((1,3,0,0), LOG_BETA[2]), "ZnCl4--": ((1,4,0,0), LOG_BETA[3]),
    "Zn(NH3)2+": ((1,0,1,0), LOG_BETA[4]), "Zn(NH3)2^2+": ((1,0,2,0), LOG_BETA[5]),
    "Zn(NH3)3^2+": ((1,0,3,0), LOG_BETA[6]), "Zn(NH3)4^2+": ((1,0,4,0), LOG_BETA[7]),
    "ZnCl3(NH3)-": ((1,3,1,0), LOG_BETA[8]), "ZnCl(NH3)3+": ((1,1,3,0), LOG_BETA[9]),
    "H+": ((0,0,0,1), 0.0), "OH-": ((0,0,0,-1), -LOG_BETA[11]),
    "NH4+": ((0,0,1,1), LOG_BETA[10]), "NH3": ((0,0,1,0), 0.0),
    "Cl-": ((0,1,0,0), 0.0), "Zn2+": ((1,0,0,0), 0.0),
}
E = np.array([_EXP[s[0]][0] for s in SPECIES], float)            # (16,4)
LNB = np.array([_EXP[s[0]][1] for s in SPECIES], float)*LN10      # (16,)

# pesos das quantidades conservadas: Zn, Cl, NH3~, NH4~
W = np.zeros((4, NS))
for nm in ["ZnCl+","ZnCl2","ZnCl3-","ZnCl4--","Zn(NH3)2+","Zn(NH3)2^2+","Zn(NH3)3^2+",
           "Zn(NH3)4^2+","ZnCl3(NH3)-","ZnCl(NH3)3+","Zn2+"]:
    W[0, IDX[nm]] = 1.0
for nm, v in [("ZnCl+",1),("ZnCl2",2),("ZnCl3-",3),("ZnCl4--",4),("ZnCl3(NH3)-",3),("ZnCl(NH3)3+",1),("Cl-",1)]:
    W[1, IDX[nm]] = v
for nm, v in [("Zn(NH3)2+",1),("Zn(NH3)2^2+",2),("Zn(NH3)3^2+",3),("Zn(NH3)4^2+",4),
              ("ZnCl3(NH3)-",1),("ZnCl(NH3)3+",3),("NH3",1),("H+",-1),("OH-",1)]:
    W[2, IDX[nm]] = v
for nm, v in [("NH4+",1),("H+",1),("OH-",-1)]:
    W[3, IDX[nm]] = v
WABS = np.abs(W)
TAU = W[[0, 2, 3], :].T.copy()            # (16,3): contribuição de cada espécie às quasipartículas
DT_RHS = np.array([[1,0,0],[2,0,1],[0,1,0],[0,0,1]], float)   # dT/dQ

# ---------------- sólidos ----------------
SOLIDS = ["ZnO", "ZnCl2_4ZnOH2_H2O", "ZnOH2", "ZnCl2_2NH3"]    # mesma ordem do final5
LOG_KSP = {"ZnO": -16.7, "ZnCl2_4ZnOH2_H2O": -14.2, "ZnOH2": -17.0, "ZnCl2_2NH3": -6.42}
# CALIBRAÇÃO: com a Tabela 2 literal, o eletrólito C (0.5 M Zn, pH 8) fica ~10x supersaturado em Zn(OH)2, mas o artigo
# afirma que C está no limite de solubilidade. Deslocar log10 Ksp dos hidróxidos em +1.0 reproduz essa afirmação e as
# demais (1o sólido em A = simonkolleite; em B = ZnCl2.2NH3). Use 0.0 para a Tabela 2 literal.
KSP_SHIFT_OH = 1.0
_SHIFT = {"ZnO": KSP_SHIFT_OH, "ZnOH2": KSP_SHIFT_OH, "ZnCl2_4ZnOH2_H2O": KSP_SHIFT_OH, "ZnCl2_2NH3": 0.0}
# SI = log10 Q - log10 Ksp, com Q por Zn: Zn*OH^m*Cl^c (simonkolleite: m=1.6, c=0.4) — ver PREMISSAS
def _si_terms():
    a, k = {}, {}
    # OH = 1/(beta_w h): ln OH = -ln beta_w - ln h
    for nm, (ez, ecl, enh, eoh) in {"ZnO": (1,0,0,2), "ZnOH2": (1,0,0,2),
                                    "ZnCl2_4ZnOH2_H2O": (1,0.4,0,1.6), "ZnCl2_2NH3": (1,2,2,0)}.items():
        a[nm] = np.array([ez, ecl, enh, -eoh], float)
        k[nm] = -eoh*LOG_BETA[11] - (LOG_KSP[nm] + _SHIFT[nm])   # constante em log10
    return a, k
SI_A, SI_K = _si_terms()
# efeito de cada sólido (por mol de sólido) em Q=[Zn~,NH3~,NH4~] — Tabela 5
NU = {"ZnO": np.array([-1,-2,2.]), "ZnOH2": np.array([-1,-2,2.]),
      "ZnCl2_4ZnOH2_H2O": np.array([-5,-8,8.]), "ZnCl2_2NH3": np.array([-1,-2,0.])}
# volumes molares (m3/mol): ZnO 14.5, Zn(OH)2 32.6, simonkolleite ~172, ZnCl2.2NH3 ~85 (estimado)
VM = {"ZnO": 14.5e-6, "ZnOH2": 32.6e-6, "ZnCl2_4ZnOH2_H2O": 172.0e-6, "ZnCl2_2NH3": 85.0e-6}
V_ZN = 9.16e-6
ZT_FLOOR = 1e-3          # mol/m3 (1 uM): piso numérico do Zn~; desvio de massa de Zn <= 1e-5 relativo

# ---------------- especiação vetorizada ----------------
def speciate(Q, y0=None, tol=1e-11, maxit=80, want_dy=False):
    """Q (n,3) mol/m3 -> y=ln([Zn2+],[Cl-],[NH3],[H+]) em mol/L, c (n,16) mol/m3, S (n,16,3) dc/dQ, ok (n,)."""
    Q = np.atleast_2d(np.asarray(Q, float))
    n = len(Q)
    q = Q/1000.0
    Zt, A3, A4 = q[:,0], q[:,1], q[:,2]
    Clt = 2*Zt + A4
    tgt = np.stack([Zt, Clt, A3, A4], 1)
    if y0 is None:
        y = _cold_start(Zt, Clt, A3, A4)
    else:
        y = np.array(y0, float).reshape(n,4).copy()

    def resid(yy):
        c = np.exp(LNB + yy@E.T)
        Fv = c@W.T - tgt
        mag = c@WABS.T + np.abs(tgt) + 1e-300
        return c, Fv, mag, np.max(np.abs(Fv)/mag, axis=1)

    c, Fv, mag, err = resid(y)
    for it in range(maxit):
        act = err > tol
        if not act.any():
            break
        J = (c[act][:,None,:]*W[None,:,:]) @ E            # (m,4,4)
        sc = mag[act]
        Jn = J/sc[:,:,None]
        rhs = -Fv[act]/sc
        try:
            dy = np.linalg.solve(Jn, rhs[:,:,None])[:,:,0]
        except np.linalg.LinAlgError:
            dy = np.linalg.lstsq(Jn[0], rhs[0], rcond=None)[0][None,:].repeat(act.sum(),0)
        dy = np.clip(dy, -3.0, 3.0)
        ia = np.where(act)[0]
        best_err = err[ia].copy(); best_y = y[ia].copy()
        done = np.zeros(len(ia), bool)
        for damp in (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03):
            yt = y[ia] + damp*dy
            ct, Ft, mt, et = resid_sub(yt, tgt[ia])
            better = (et < best_err) & ~done
            best_err[better] = et[better]; best_y[better] = yt[better]; done |= better
            if done.all():
                break
        y[ia] = best_y
        c, Fv, mag, err = resid(y)
    ok = err < 1e-8
    S = None
    dydq = None
    if want_dy or True:
        J = (c[:,None,:]*W[None,:,:]) @ E
        try:
            dydq = np.linalg.solve(J, np.broadcast_to(DT_RHS, (n,4,3)).copy())
        except np.linalg.LinAlgError:
            dydq = np.zeros((n,4,3)); ok &= False
        S = c[:,:,None]*np.einsum('sk,nkp->nsp', E, dydq)
    return y, c*1000.0, S, ok, dydq

def _cold_start(Zt, Clt, A3, A4):
    n = len(Zt)
    nh3 = np.maximum(A3 - 3.5*Zt, 0.02*np.abs(A3) + 1e-7)
    h = np.clip(np.maximum(A4, 1e-9)/(10**LOG_BETA[10]*nh3), 1e-14, 1.0)
    y = np.empty((n,4))
    y[:,0] = np.log(np.maximum(Zt*0.05, 1e-14))
    y[:,1] = np.log(np.maximum(Clt*0.6, 1e-14))
    y[:,2] = np.log(nh3)
    y[:,3] = np.log(h)
    return y

def resid_sub(yy, tgt):
    c = np.exp(LNB + yy@E.T)
    Fv = c@W.T - tgt
    mag = c@WABS.T + np.abs(tgt) + 1e-300
    return c, Fv, mag, np.max(np.abs(Fv)/mag, axis=1)

def si_all(y):
    """SI (n,4) na ordem SOLIDS, a partir de y."""
    out = np.empty((len(y), 4))
    for k, nm in enumerate(SOLIDS):
        out[:,k] = (y@SI_A[nm])/LN10 + SI_K[nm]
    return out

def init_composition(zncl2_M, nh4cl_M, pH):
    """Q0 (mol/m3) com NH4OH ajustado ao pH: Zn~, NH4~ fixos; NH3~ resulta do pH."""
    Zt = zncl2_M; A4 = nh4cl_M; Clt = 2*Zt + A4; h = 10.0**(-pH)
    y = np.array([math.log(Zt*0.05), math.log(Clt*0.6), math.log(max(A4,1e-6)/(10**LOG_BETA[10]*h))])
    rows = [0, 1, 3]
    for it in range(200):
        yy = np.array([y[0], y[1], y[2], math.log(h)])
        c = np.exp(LNB + E@yy)
        Fv = np.array([W[0]@c - Zt, W[1]@c - Clt, W[3]@c - A4])
        mag = np.array([WABS[0]@c + Zt, WABS[1]@c + Clt, WABS[3]@c + A4])
        if np.max(np.abs(Fv)/mag) < 1e-12:
            break
        J = (c[None,:]*W[rows])@E[:, :3]
        dy = np.linalg.solve(J/mag[:,None], -Fv/mag)
        y = y + np.clip(dy, -2, 2)
    else:
        raise ArithmeticError("init_composition não convergiu")
    A3 = float(W[2]@c)
    return np.array([Zt, A3, A4])*1000.0

# ---------------- modelo 1D ----------------
class Model:
    def __init__(self, setup=None):
        s = dict(setup or {})
        self.setup = s
        mm = 1e-3
        self.zn_model = s.get("zn_model", "paste")           # "paste" | "foil"
        paste = self.zn_model == "paste"
        self.L_ZN = float(s.get("zn_mm", 5.0))*mm if paste else 0.0
        self.L_SEP = float(s.get("separator_mm", 0.5))*mm
        self.L_GDE = float(s.get("gde_mm", 0.5))*mm
        n_zn = int(s.get("n_zn", 16)) if paste else 0
        n_sep = int(s.get("n_sep", 8 if self.L_SEP < 2*mm else 24))
        n_gde = int(s.get("n_gde", 4))
        # malha: Zn uniforme, separador com refino nas duas pontas, GDE uniforme
        ss = np.linspace(0, 1, n_sep+1)
        g = 1.6
        e_sep = self.L_SEP*(1 + np.tanh(g*(2*ss-1))/math.tanh(g))/2
        dx = []
        if n_zn: dx += [self.L_ZN/n_zn]*n_zn
        dx += list(np.diff(e_sep))
        dx += [self.L_GDE/n_gde]*n_gde
        self.dx = np.array(dx)
        self.region = ["ZN"]*n_zn + ["SEP"]*n_sep + ["GDE"]*n_gde
        self.N = len(self.dx)
        self.n_zn, self.n_sep, self.n_gde = n_zn, n_sep, n_gde
        self.iz = np.arange(0, n_zn)
        self.ig = np.arange(n_zn+n_sep, self.N)
        self.dxf = 0.5*(self.dx[:-1] + self.dx[1:])
        self.x = np.cumsum(self.dx) - 0.5*self.dx
        # porosidade do eletrólito (estática) e Bruggeman
        self.beta = float(s.get("bruggeman", 1.5))
        eps = np.empty(self.N)
        for i, r in enumerate(self.region):
            eps[i] = {"ZN": float(s.get("eps_e_zn", 0.30)), "SEP": float(s.get("eps_e_sep", 0.90)),
                      "GDE": float(s.get("eps_e_gde", 0.25))}[r]
        self.eps_e = eps
        self.epsb_f = 0.5*(eps[:-1]**self.beta + eps[1:]**self.beta)
        # zinco
        self.eps_zn0 = float(s.get("eps_zn0", 0.40))
        self.d_zn = float(s.get("d_zn_um", 100.0))*1e-6
        self.delta_c = float(s.get("delta_pass_um", 5.0))*1e-6
        self.eps_zn = np.zeros(self.N)
        if paste:
            self.eps_zn[self.iz] = self.eps_zn0
        self.foil_rf = float(s.get("foil_rf", 1.0))
        # cinética (calibrável) — ver PREMISSAS
        self.i0_zn = float(s.get("i0_zn", 2*F*8.0e-6))        # A/m2 (k0 do SI x 2F)
        self.i0_orr = float(s.get("i0_orr", 1.86e-6))         # A/m2 de área ativa (calibrado: V=0.90 V @1 mA/cm2, eletrólito C)
        self.i0_oer = float(s.get("i0_oer", 4.49e-3))         # A/m2 de área ativa (calibrado: V=2.25 V @2.5 mA/cm2, célula A)
        self.alpha_zn = 0.5
        self.alpha_orr = float(s.get("alpha_orr", 0.5))
        self.alpha_oer = float(s.get("alpha_oer", 0.5))
        self.A_gde = float(s.get("A_gde_sp", 4.5e3))
        self.pO2 = float(s.get("pO2_atm", 0.21))
        self.solids_on = [k for k in SOLIDS if k in s.get("solids", SOLIDS)]
        # estado
        zt = float(s.get("zncl2_M", 0.5)); nh = float(s.get("nh4cl_M", 1.6))
        self.initial_pH = float(s.get("initial_pH", 8.0))
        Q0 = init_composition(zt, nh, self.initial_pH)
        self.Q = np.tile(Q0, (self.N, 1))
        self.s = np.zeros((4, self.N))                    # sólidos, mol/m3 de volume total
        self.y = None
        self._respec()
        self.t = 0.0; self.cap = 0.0; self.energy = 0.0   # Ah/m2, Wh/m2
        self.V = float('nan')
        self.dt = float(s.get("dt0", 1.0)); self.dt_max = float(s.get("dt_max", 120.0))
        self.dpH_max = float(s.get("dpH_max", 1.0))
        self.history = []
        self.n_steps = 0; self.n_rej = 0
        self.stop_flag = False
        _g = float(self.pH()[self.ig].mean())
        self.pH_min_run = _g; self.pH_max_run = _g
        self._equilibrate_initial(zt, nh)
        self.inv0 = self.invariants()

    # ----- utilitários -----
    def _respec(self):
        y, c, S, ok, dydq = speciate(self.Q, self.y)
        if not ok.all():
            y2, c2, S2, ok2, dydq2 = speciate(self.Q[~ok], None)
            y[~ok], c[~ok], S[~ok], ok[~ok] = y2, c2, S2, ok2
            if not ok.all():
                raise ArithmeticError("especiação não convergiu")
        self.y, self.c, self.S = y, c, S

    def pH(self):
        return -np.log10(np.maximum(self.c[:, IH]/1000.0, 1e-30))

    def solids_vol(self):
        return sum(VM[k]*self.s[i] for i, k in enumerate(SOLIDS))

    def invariants(self):
        dx = self.dx; ee = self.eps_e
        zn = np.dot(dx, ee*self.Q[:,0]) + np.dot(dx, self.eps_zn)/V_ZN
        cl = np.dot(dx, ee*(2*self.Q[:,0] + self.Q[:,2]))
        nn = np.dot(dx, ee*(self.Q[:,1] + self.Q[:,2]))
        s = {k: np.dot(dx, self.s[i]) for i, k in enumerate(SOLIDS)}
        zn += 5*s["ZnCl2_4ZnOH2_H2O"] + s["ZnO"] + s["ZnOH2"] + s["ZnCl2_2NH3"]
        cl += 2*s["ZnCl2_4ZnOH2_H2O"] + 2*s["ZnCl2_2NH3"]
        nn += 2*s["ZnCl2_2NH3"]
        return np.array([zn, cl, nn])

    def _equilibrate_initial(self, zt, nh):
        """Eletrólito preparado: precipita o excedente e reajusta o pH com NH4OH (só NH3~ muda)."""
        for it in range(12):
            self._flash(self.Q, self.s); self._respec()
            q0 = self.Q[0]/1000.0
            A3 = init_composition(q0[0], q0[2], self.initial_pH)[1]
            dA3 = A3 - self.Q[0,1]
            self.Q[:,1] += dA3
            self._respec()
            if abs(dA3) < 1e-6*max(1.0, abs(A3)) and np.max(np.abs(self.pH()-self.initial_pH)) < 1e-6:
                break
        self._flash(self.Q, self.s); self._respec()

    # ----- precipitação / dissolução em equilíbrio -----
    def _flash(self, Q, s, tol_si=2e-6):
        ee = self.eps_e
        y = self.y if self.y is not None else speciate(Q)[0]
        SIc = si_all(y)
        for sweep in range(3):
            moved = False
            for k, nm in enumerate(SOLIDS):
                if nm not in self.solids_on:
                    continue
                SI = SIc[:, k]
                have = s[k] > 1e-9
                idx = np.where((SI > tol_si) | ((SI < -tol_si) & have))[0]
                if len(idx) == 0:
                    continue
                nu = NU[nm]
                Qi = Q[idx].copy(); yi = y[idx].copy()
                se = s[k][idx]/ee[idx]                           # sólido na base eletrólito
                hi = np.where(SI[idx] > 0, 0.99*Qi[:,0]/abs(nu[0]), 0.0)
                lo = np.where(SI[idx] < 0, -se, 0.0)
                xi = np.zeros(len(idx))
                a = SI_A[nm]
                for it in range(14):
                    Qt = Qi + xi[:,None]*nu
                    yt, ct, St, okt, dydq = speciate(Qt, yi)
                    f = (yt@a)/LN10 + SI_K[nm]
                    f = np.where(okt, f, -1.0)                  # falha -> trata como supersaturação negativa
                    lo = np.where(f > 0, xi, lo); hi = np.where(f <= 0, xi, hi)
                    d = np.einsum('k,nkp,p->n', a, dydq, nu)/LN10/1000.0
                    with np.errstate(all='ignore'):
                        new = xi - f/d
                    bad = ~np.isfinite(new) | (new <= lo) | (new >= hi)
                    new = np.where(bad, 0.5*(lo+hi), new)
                    conv = (np.abs(f) < 2e-7) | ((hi - lo) < 1e-12*np.maximum(1.0, np.abs(hi)+np.abs(lo)))
                    xi = np.where(conv, xi, new)
                    yi = np.where(okt[:,None], yt, yi)
                    if conv.all():
                        break
                Q[idx] = Qi + xi[:,None]*nu
                yi, _c, _S, _ok, _d = speciate(Q[idx], yi)
                s[k][idx] = np.maximum(s[k][idx] + ee[idx]*xi, 0.0)
                moved = True
                y[idx] = yi
                SIc[idx] = si_all(yi)
            if not moved:
                break
        self.y = y
        self._flash_moved = moved

    # ----- zinco: área ativa e passivação -----
    def _zn_area(self, I):
        n = self.N; az = np.zeros(n)
        if self.zn_model == "foil":
            az[self.n_zn] = self.foil_rf/self.dx[self.n_zn]
            return az, np.ones(n)
        iz = self.iz
        eps = np.maximum(self.eps_zn[iz], 1e-3)
        A = 6.0/self.d_zn*self.eps_zn0**(1/3)*eps**(2/3)
        npart = 6.0*self.eps_zn0/(math.pi*self.d_zn**3)
        Vs = self.solids_vol()[iz]
        r = 0.5*self.d_zn*(eps/self.eps_zn0)**(1/3)
        Rr = (r**3 + 3.0*Vs/(4*math.pi*npart))**(1/3)
        delta = Rr - r
        f = np.clip(1.0 - (delta/self.delta_c)**2, 0.0, 1.0) if I > 0 else np.ones(len(iz))
        az[iz] = A*f
        self._delta = delta
        self._fmean = float(np.mean(f))
        return az, f

    # ----- potencial do eletrólito (inclui potencial de difusão) -----
    def _phi_e(self, c, jf):
        cf = 0.5*(c[:-1] + c[1:])
        gr = (c[1:] - c[:-1])/self.dxf[:,None]
        G = cf@Z2D
        dif = gr@ZD
        dphi = -RTF*((jf[1:-1]/F)/self.epsb_f + dif)/np.maximum(G, 1e-12)*self.dxf
        return np.concatenate([[0.0], np.cumsum(dphi)])

    def _solve_phis(self, base, w, ai, ac, target):
        """Σ w (exp(ai(φ-b)/RTF) - exp(-ac(φ-b)/RTF)) = target ; Newton monótono."""
        m = w > 0
        if not m.any():
            return None
        b = base[m]; ww = w[m]
        phi = float(np.dot(ww, b)/ww.sum())
        for it in range(200):
            e1 = np.exp(np.clip(ai*(phi-b)/RTF, -80, 80)); e2 = np.exp(np.clip(-ac*(phi-b)/RTF, -80, 80))
            g = np.dot(ww, e1 - e2) - target
            dg = np.dot(ww, ai/RTF*e1 + ac/RTF*e2)
            if abs(g) < 1e-12*max(abs(target), 1e-12):
                break
            step = np.clip(-g/dg, -0.25, 0.25)
            phi += step
        return phi

    def kinetics(self, I):
        """I em A/m2 (>0 descarga). Devolve rz, rg (A/m3), V, diagnósticos."""
        c = self.c; cM = c/1000.0; dx = self.dx; n = self.N
        Uz = -0.762 + RTF/2*np.log(np.maximum(cM[:,IZN], 1e-30))
        Ug = 1.229 + RTF/4*math.log(self.pO2) + RTF*np.log(np.maximum(cM[:,IH], 1e-30))
        az, fp = self._zn_area(I)
        ag = np.zeros(n); ag[self.ig] = self.A_gde
        zmask = az > 0
        if I == 0.0:
            wz = az*dx; wg = ag*dx
            V = float(np.dot(wg, Ug)/wg.sum() - np.dot(wz, Uz)/wz.sum())
            return np.zeros(n), np.zeros(n), V, {"eta_z": 0.0, "eta_g": 0.0}
        if az.sum() <= 0:
            return np.zeros(n), np.zeros(n), -1e9, {"passivated": True}
        i0g = self.i0_orr if I > 0 else self.i0_oer
        wz0 = az*dx; wg0 = ag*dx
        rz = I*wz0/wz0.sum(); rg = -I*wg0/wg0.sum()
        for it in range(3):
            jf = np.concatenate([[0.0], np.cumsum((rz + rg)*dx)])
            phi = self._phi_e(c, jf)
            wz = self.i0_zn*az*dx; wg = i0g*ag*dx
            ag_ = self.alpha_orr if I > 0 else self.alpha_oer
            pz = self._solve_phis(phi + Uz, wz, self.alpha_zn, self.alpha_zn, I)
            pg = self._solve_phis(phi + Ug, wg, ag_, ag_, -I)
            if pz is None or pg is None:
                return np.zeros(n), np.zeros(n), -1e9, {"passivated": True}
            ez = np.clip((pz - phi - Uz)/RTF*self.alpha_zn, -80, 80)
            eg = np.clip((pg - phi - Ug)/RTF*ag_, -80, 80)
            rz = self.i0_zn*az*(np.exp(ez) - np.exp(-ez))
            rg = i0g*ag*(np.exp(eg) - np.exp(-eg))
        V = pg - pz
        return rz, rg, float(V), {"eta_z": float(np.dot(az*dx, pz-phi-Uz)/max((az*dx).sum(),1e-30)),
                                  "eta_g": float(np.dot(ag*dx, pg-phi-Ug)/(ag*dx).sum()),
                                  "ohm": float(phi[self.ig].mean() - phi[:max(self.n_zn,1)].mean()) if I != 0 else 0.0}

    # ----- um passo -----
    def _transport_matrix(self, jf):
        c, S = self.c, self.S
        n = self.N
        cf = 0.5*(c[:-1] + c[1:]); Sf = 0.5*(S[:-1] + S[1:])
        G = np.maximum(cf@Z2D, 1e-12)
        A = np.einsum('iq,i,fip->fqp', TAU, D_SP, Sf)
        Bv = np.einsum('i,fip->fp', ZD, Sf)
        M = np.einsum('iq,fi->fq', TAU*ZD[:,None], cf)/G[:,None]
        K = self.epsb_f[:,None,None]*(A - M[:,:,None]*Bv[:,None,:])/self.dxf[:,None,None]
        m = M*(jf[1:-1]/F)[:,None]
        return K, m

    def advance(self, I, dt, kin=None):
        n = self.N; dx = self.dx; ee = self.eps_e
        rz, rg, V, info = kin if kin is not None else self.kinetics(I)
        src = np.zeros((n, 3))
        src[:,0] = rz/(2*F)
        src[:,1] = -rg/F
        src[:,2] = rg/F
        jf = np.concatenate([[0.0], np.cumsum((rz + rg)*dx)])
        K, m = self._transport_matrix(jf)
        Am = np.zeros((3*n, 3*n)); rhs = np.zeros(3*n)
        Id = np.eye(3)
        for i in range(n):
            sl = slice(3*i, 3*i+3)
            Dg = ee[i]*dx[i]/dt*Id
            fl = np.zeros(3); fr = np.zeros(3)
            if i > 0:
                Dg = Dg + K[i-1]; Am[sl, 3*i-3:3*i] = -K[i-1]; fl = m[i-1]
            if i < n-1:
                Dg = Dg + K[i];   Am[sl, 3*i+3:3*i+6] = -K[i];   fr = m[i]
            Am[sl, sl] = Dg
            rhs[sl] = ee[i]*dx[i]/dt*self.Q[i] + (fl - fr) + dx[i]*src[i]
        Qn = np.linalg.solve(Am, rhs).reshape(n, 3)
        return Qn, V, info, rz

    def step(self, I):
        """Avança um passo adaptativo. Devolve (V, info)."""
        kin = self.kinetics(I)
        V = kin[2]
        if V < -1e8:
            return V, kin[3]
        while True:
            dt = self.dt
            snap = (self.Q.copy(), self.s.copy(), self.eps_zn.copy(), self.y.copy())
            ph0 = self.pH()
            try:
                Qn, V, info, rz = self.advance(I, dt, kin)
                if not np.all(np.isfinite(Qn)):
                    raise ArithmeticError("Q não finito")
                Qn[:,0] = np.maximum(Qn[:,0], ZT_FLOOR)       # piso numérico de Zn~ (volume esgotado)
                self.Q = Qn
                if self.zn_model == "paste":
                    self.eps_zn[self.iz] = np.maximum(self.eps_zn[self.iz] - V_ZN*rz[self.iz]/(2*F)*dt, 0.0)
                self.y = np.clip(self.y, -60, 20)
                y, c, S, ok, _ = speciate(self.Q, self.y)
                if not ok.all():
                    y2, c2, S2, ok2, _ = speciate(self.Q[~ok], None)
                    if not ok2.all():
                        raise ArithmeticError("especiação")
                    y[~ok] = y2
                self.y = y
                self._flash(self.Q, self.s)
                self._respec()
                dph = float(np.max(np.abs(self.pH() - ph0)))
                if dph > self.dpH_max and dt > 1e-3:
                    raise ArithmeticError("dpH")
            except (ArithmeticError, np.linalg.LinAlgError):
                self.Q, self.s, self.eps_zn, self.y = snap
                self._respec()
                self.dt = max(dt*0.5, 1e-4); self.n_rej += 1
                if dt <= 1e-4:
                    raise
                continue
            break
        self.t += dt; self.n_steps += 1
        self.cap += I*dt/3600.0
        self.energy += I*V*dt/3600.0
        self.V = V
        _g = float(self.pH()[self.ig].mean())
        self.pH_min_run = min(self.pH_min_run, _g); self.pH_max_run = max(self.pH_max_run, _g)
        self.dt = min(self.dt_max, dt*(1.4 if dph < 0.5*self.dpH_max else 1.0))
        return V, info

    # ----- protocolo -----
    def run_protocol(self, segments, rec_every_s=60.0, callback=None, max_steps=2_000_000, cb_wall_s=0.25):
        """segments: lista de dict(mode='discharge'|'charge'|'rest', I_mA_cm2, t_max_h, V_cut, q_max_mAh_cm2).
        callback(model, k, nseg, frac, V, cpu_s) é chamado a cada cb_wall_s s de relógio; stop_flag interrompe."""
        t0 = time.time(); next_rec = self.t; last_cb = 0.0; self.stop_flag = False
        out = []
        for seg in segments:
            mode = seg.get("mode", "discharge")
            I = seg.get("I_mA_cm2", 1.0)*10.0*(1 if mode == "discharge" else -1 if mode == "charge" else 0)
            t_end = self.t + seg.get("t_max_h", 1e9)*3600.0
            vcut = seg.get("V_cut", 0.6 if mode == "discharge" else 2.6 if mode == "charge" else None)
            q_end = self.cap + (seg.get("q_max_mAh_cm2", 1e9)*10.0)*(1 if I >= 0 else -1)   # mAh/cm2 -> Ah/m2: x10
            reason = "tempo"
            seg_t0 = self.t; seg_c0 = self.cap; seg_e0 = self.energy
            while True:
                if self.stop_flag: reason = "interrompido"; break
                if self.t >= t_end: reason = "tempo"; break
                if I > 0 and self.cap >= q_end: reason = "capacidade"; break
                if I < 0 and self.cap <= q_end: reason = "capacidade"; break
                dt_left = t_end - self.t
                self.dt = min(self.dt, dt_left)
                V, info = self.step(I)
                if V < -1e8:
                    reason = "passivação (Zn sem área ativa)"; break
                if vcut is not None and ((I > 0 and V < vcut) or (I < 0 and V > vcut)):
                    reason = f"V_cut {vcut} V"; break
                if self.t >= next_rec:
                    self.history.append(self.record(I, V, info)); next_rec = self.t + rec_every_s
                if callback and time.time() - last_cb >= cb_wall_s:
                    last_cb = time.time()
                    frac = min(1.0, (self.t - seg_t0)/max(t_end - seg_t0, 1e-9))
                    callback(self, len(out), len(segments), frac, V, time.time() - t0)
                if self.n_steps > max_steps: reason = "max_steps"; break
            if self.stop_flag: out.append({"mode": mode, "I_mA_cm2": I/10.0, "t_h": (self.t-seg_t0)/3600.0, "cap_mAh_cm2": (self.cap-seg_c0)/10.0, "E_Wh_m2": self.energy-seg_e0, "V_end": self.V, "reason": reason, "pH_gde": float(self.pH()[self.ig].mean()), "pH_zn": float(self.pH()[:self.n_zn+1].mean())}); break
            out.append({"mode": mode, "I_mA_cm2": I/10.0, "t_h": (self.t-seg_t0)/3600.0,
                        "cap_mAh_cm2": (self.cap-seg_c0)/10.0, "E_Wh_m2": self.energy-seg_e0,
                        "V_end": self.V, "reason": reason, "pH_gde": float(self.pH()[self.ig].mean()),
                        "pH_zn": float(self.pH()[:self.n_zn+1].mean())})
        self.cpu_s = time.time() - t0
        return out

    def record(self, I, V, info):
        ph = self.pH()
        return {"t_h": self.t/3600.0, "V": V, "I": I/10.0, "pH_gde": float(ph[self.ig].mean()),
                "pH_min": float(ph.min()), "pH_max": float(ph.max()),
                "pH_zn": float(ph[:max(self.n_zn,1)].mean()),
                "solids_mol_m2": [float(np.dot(self.dx, self.s[k])) for k in range(4)],
                "eps_zn_mean": float(self.eps_zn[self.iz].mean()) if self.n_zn else 1.0,
                "cap_mAh_cm2": self.cap/10.0,
                "delta_um": (float(np.max(self._delta))*1e6 if (self.n_zn and hasattr(self, "_delta")) else 0.0),
                "f_pass": float(getattr(self, "_fmean", 1.0)),
                "dt": self.dt, **{k: v for k, v in info.items() if isinstance(v, float)}}

    def snapshot(self):
        """Cópia leve do estado espacial para gráficos (chamar na mesma thread que avança o modelo)."""
        zf = lambda idx: {SPECIES[k][0]: float(v) for k, v in enumerate(
            (self.c[idx]*W[0]).mean(0)/max((self.c[idx]*W[0]).mean(0).sum(), 1e-30))}
        return {"x": list(self.x*1000.0), "edges": list(np.cumsum(self.dx)*1000.0),
                "n_zn": self.n_zn, "n_sep": self.n_sep, "pH": [float(v) for v in self.pH()],
                "Zt": [float(v) for v in self.Q[:,0]/1000.0], "A3": [float(v) for v in self.Q[:,1]/1000.0],
                "A4": [float(v) for v in self.Q[:,2]/1000.0], "z_gde": zf(self.ig),
                "z_zn": zf(np.arange(0, max(self.n_zn, 1)))}

    def energy_density_WhL(self):
        """Wh por litro de volume total da célula (Zn + separador + GDE)."""
        return self.energy/self.dx.sum()/1000.0

# PREMISSAS (estáticas salvo indicação):
#  - Concentrações no lugar de atividades (I ~ 3-6 M): erro conhecido nas constantes e no Nernst.
#  - log10 Ksp dos hidróxidos deslocado em +1.0 (KSP_SHIFT_OH): calibração, não valor da Tabela 2.
#  - Ksp simonkolleite em base por-Zn: Zn*OH^1.6*Cl^0.4 (reproduz o cruzamento pH~7.5 a Cl~3 M do artigo).
#  - Volume molar ZnCl2.2NH3 ~85 cm3/mol (estimado); eps_e fixo por região; sem inundação do GDE e sem balanço de H2O.
#  - Passivação: casca de sólido sobre partícula esférica de Zn; ativa=1-(delta/5um)^2 (forma assumida).
#  - O2 no GDE: pO2 constante (dinâmico apenas se o usuário alterar pO2_atm).
#  - i0_orr / i0_oer: calibrados em curvas digitalizadas do artigo (Fig. 5/8), não medidos.

# ===== FIM DO NÚCLEO MATEMÁTICO =====
# ============================================================
# VERBETES TÉCNICOS — apenas grandezas e regras efetivamente presentes no código v2
# ============================================================
S01 = "01 — ARQUITETURA E NUMÉRICA"
S02 = "02 — CONSTANTES, CONDIÇÕES E GEOMETRIA"
S03 = "03 — ESPÉCIES QUÍMICAS"
S04 = "04 — QUASIPARTÍCULAS E ESPECIAÇÃO"
S05 = "05 — POTENCIAIS E CINÉTICA"
S06 = "06 — TRANSPORTE"
S07 = "07 — SÓLIDOS E PRECIPITAÇÃO"
S08 = "08 — ZINCO, PASSIVAÇÃO E GDE"
S09 = "09 — PROTOCOLO E RESULTADOS"
S10 = "10 — VARREDURA DE ESPESSURA DE Zn"
S11 = "11 — PREMISSAS E CALIBRAÇÕES"
S12 = "12 — REFERÊNCIAS BIBLIOGRÁFICAS"

GLOSSARY_TERMS = [
# ---------------- 01 ----------------
(S01, "Modelo 1D por volumes finitos (v2)", "Domínio unidimensional Zn | separador | GDE dividido em volumes finitos. O estado de cada volume são as três quasipartículas Zn~, NH3~ e NH4~; as 16 espécies são reconstruídas por especiação a cada passo."),
(S01, "Domínio físico e regiões", "Vetor region com os rótulos ZN, SEP e GDE, um por volume. Com Zn em folha (foil) não existem volumes ZN: o separador começa na superfície do metal."),
(S01, "Malha", "Zn uniforme (16 volumes, ou max(16, 4×mm) na varredura), separador com refino em tanh nas duas pontas (8 volumes se < 2 mm, senão 24) e GDE uniforme (4 volumes). dx guarda a espessura de cada volume em metros."),
(S01, "Volume finito", "Unidade espacial da malha. Cada volume guarda Q, os sólidos, a fração de Zn (eps_zn) e a fração de eletrólito (eps_e)."),
(S01, "Vetor de estado Q", "Matriz (N,3) em mol/m³ de eletrólito com colunas [Zn~, NH3~, NH4~]. É tudo o que o transporte e as fontes atualizam; Cl e N totais são derivados."),
(S01, "Passo de tempo adaptativo", "Não há 'passos × dt'. dt começa em 1 s, cresce ×1,4 quando a variação de pH do passo é menor que metade de dpH_max (1,0) e é limitado por dt_max (120 s; 300 s na varredura)."),
(S01, "Passo rejeitado", "Passo desfeito e repetido com dt/2 quando a especiação não converge, aparece valor não finito ou o pH de algum volume muda mais que dpH_max. O estado anterior é restaurado. Com dt ≤ 1e-4 s o erro é propagado."),
(S01, "Passos aceitos / rejeitados", "Contadores do relatório: n_steps e n_rej. Muitas rejeições indicam região rígida (frente ácida, esgotamento de Zn em um volume)."),
(S01, "Esquema semi-implícito", "A difusão é implícita (Euler regressivo) com coeficientes congelados no início do passo; fontes e migração são explícitas. Não há limite de estabilidade do tipo dx²/2D."),
(S01, "Matriz bloco-tridiagonal", "Sistema 3N×3N (3 quasipartículas por volume) montado a cada passo e resolvido por numpy.linalg.solve. Conserva Zn~, NH3~ e NH4~ por construção (fluxos se anulam entre volumes vizinhos)."),
(S01, "Divisão de operadores", "Ordem de um passo: cinética (corrente e V) → transporte implícito com fontes → atualização de eps_zn → especiação → precipitação/dissolução → especiação final. Erro de divisão controlado pelo limite de variação de pH."),
(S01, "Tempo de simulação", "Tempo físico representado (model.t, em s ou h). É a soma dos dt aceitos."),
(S01, "Tempo computacional", "Tempo de CPU da execução (cpu_s). Medida de desempenho, não variável física."),
(S01, "Estado final", "Valores armazenados depois do último passo: Q, sólidos, eps_zn, pH e perfis. O gráfico usa o histórico (history), gravado a cada rec_every_s segundos simulados."),
(S01, "Invariantes e deriva", "Totais de Zn (inclui metal e sólidos), Cl e N integrados no domínio. A deriva relativa em relação ao início é impressa no relatório; deve ficar na ordem de 1e-12 (até ~1e-5 se o piso de Zn~ atuar)."),
(S01, "Piso numérico de Zn~", "ZT_FLOOR = 1e-3 mol/m³ (1 µM). Evita log(0) na especiação quando um volume esgota o Zn dissolvido, como no GDE durante a carga."),
(S01, "Interrupção (PARAR)", "stop_flag do modelo, verificada a cada passo. O segmento em curso é encerrado com motivo 'interrompido' e o relatório usa o estado atingido."),
# ---------------- 02 ----------------
(S02, "F — constante de Faraday", "96485,33212 C/mol. Converte corrente em mol de elétrons."),
(S02, "R — constante dos gases", "8,314462618 J/(mol·K)."),
(S02, "T — temperatura", "298,15 K, fixa no código (estática)."),
(S02, "RT/F", "Escala térmica do potencial, ≈ 25,7 mV. Aparece nos termos de Nernst, de Butler-Volmer e no potencial do eletrólito."),
(S02, "Condição galvanostática", "A corrente total é imposta; a tensão da célula resulta. Sinal: corrente > 0 é descarga, < 0 é carga, 0 é repouso."),
(S02, "Corrente em mA/cm² e A/m²", "A interface usa mA/cm²; o núcleo usa A/m². 1 mA/cm² = 10 A/m². Capacidade: 1 mAh/cm² = 10 Ah/m²."),
(S02, "Composição inicial", "Concentrações de ZnCl2 e NH4Cl (mol/L). Cl total = 2·Zn + NH4Cl e NH4~ = NH4Cl; NH3~ é determinado pelo pH inicial."),
(S02, "pH inicial", "Ajustado por NH4OH, que só altera NH3~ (NH4+ e OH- criados em igual quantidade não mudam NH4~)."),
(S02, "Equilíbrio inicial", "Antes de simular, o excedente supersaturado precipita e o pH é reajustado com NH3~ (até 12 iterações). Por isso o estado em t = 0 pode ter pouco sólido e Zn dissolvido ligeiramente abaixo do nominal."),
(S02, "Eletrólitos A, B, C e D", "A: 0,51 M ZnCl2 + 2,34 M NH4Cl, pH 6. B: 0,26 M + 5 M, pH 7. C: 0,5 M + 1,6 M, pH 8. D: 0,5 M + 2 M, pH 7 com só ZnO habilitado (caso idealizado do artigo)."),
(S02, "Espessuras Zn, separador e GDE", "Zn: 5 mm no padrão (só vale em 'paste'). Separador: 0,5 mm no padrão; a Célula A do artigo usa 30 mm. GDE: 0,5 mm."),
(S02, "Modelo de Zn: paste e foil", "paste: pasta de partículas esféricas de 100 µm, com fração de Zn 0,40, ocupando volumes ZN. foil: folha polida; a reação ocorre na superfície, no primeiro volume do separador, com área = fator de rugosidade / dx."),
(S02, "Fração de eletrólito eps_e", "Estática por região: Zn 0,30, separador 0,90, GDE 0,25. Multiplica o armazenamento e, elevada a 1,5, a difusão."),
(S02, "Coeficiente de Bruggeman", "Expoente 1,5 aplicado a eps_e para obter a difusão efetiva em meio poroso."),
(S02, "pO2 — pressão de O2", "0,21 atm (ar), constante. Só entra no potencial de equilíbrio do GDE; não há transporte de O2 no modelo."),
(S02, "Sólidos habilitados", "'todos' liga os 4 sólidos; 'só ZnO' reproduz o eletrólito D idealizado do artigo (ZnO como único produto)."),
# ---------------- 03 ----------------
(S03, "Espécies químicas (16)", "Zn2+, 4 complexos Zn-Cl, 4 complexos Zn-NH3, 2 ternários, H+, OH-, NH4+, NH3 e Cl-. Concentrações em mol/m³ no estado e em mol/L dentro da especiação. Ao contrário do final5, O2 e H2O não são espécies."),
(S03, "Zn2+", "Íon livre. Sua concentração é muito pequena em eletrólito amoniacal (complexos dominam) e define o potencial de equilíbrio do Zn."),
(S03, "ZnCl+, ZnCl2, ZnCl3-, ZnCl4--", "Complexos cloreto, favorecidos em pH ácido, onde NH3 livre é baixo."),
(S03, "Zn(NH3)2+, Zn(NH3)2^2+, Zn(NH3)3^2+, Zn(NH3)4^2+", "Complexos amínicos com 1, 2, 3 e 4 NH3 (nomenclatura do código). Zn(NH3)4^2+ domina em pH alto e no GDE ao fim da descarga."),
(S03, "ZnCl3(NH3)- e ZnCl(NH3)3+", "Complexos ternários Zn-Cl-NH3, importantes na faixa quase neutra."),
(S03, "H+ e OH-", "Base do pH. OH- = 1/(β_w·[H+]) com log10 β_w = 13,96."),
(S03, "NH4+ e NH3", "Par ácido-base com log10 β = 9,80 (NH3 + H+ ⇌ NH4+). NH3 também forma os complexos amínicos."),
(S03, "Cl-", "Cloreto livre. Cl total = 2·Zn~ + NH4~."),
(S03, "Complexação Zn-Cl, Zn-NH3 e ternária", "Equilíbrios de massa-ação da Tabela 1 do artigo: [complexo] = β·[Zn2+]^a·[Cl-]^b·[NH3]^c."),
(S03, "Constantes de formação β e log10 β", "12 valores de LOG_BETA (Tabela 1). Internamente usa-se ln β somado ao produto de expoentes pelas variáveis logarítmicas."),
(S03, "Concentração em mol/m³ e mol/L", "O estado e o transporte usam mol/m³; a especiação trabalha em mol/L (÷1000). As sensibilidades S são adimensionais."),
(S03, "Capacidade tamponante", "Os complexos Zn-NH3 retiram NH3 livre, permitindo que NH4+ ⇌ NH3 + H+ funcione como fonte de prótons sem deixar o pH alcalino. Entre pH ~6 e ~9,8 o pH é estável enquanto houver NH3 disponível."),
# ---------------- 04 ----------------
(S04, "Quasipartícula", "Quantidade conservada pelas reações homogêneas. Três bastam para 16 espécies e 12 reações: Zn~, NH3~ e NH4~ (Eqs. 21 a 23 do artigo)."),
(S04, "Zn~", "Zn dissolvido em todas as formas (Zn2+ + complexos Zn-Cl + Zn-NH3 + ternários). Fonte: dissolução/deposição de Zn."),
(S04, "NH3~", "[NH3] + Σ n·[Zn(NH3)n] + [ZnCl3(NH3)-] + 3·[ZnCl(NH3)3+] − [H+] + [OH-]. Fonte: ORR/OER no GDE e reações de precipitação."),
(S04, "NH4~", "[NH4+] + [H+] − [OH-]. Fonte: ORR/OER no GDE e reações de precipitação."),
(S04, "N total e Cl total derivados", "N = NH3~ + NH4~ e Cl = 2·Zn~ + NH4~ (eletroneutralidade). O Cl não é variável independente; a identidade se mantém porque a corrente do eletrólito fecha o balanço de carga."),
(S04, "Especiação", "Dado Q, resolve [Zn2+], [Cl-], [NH3] e [H+]. Newton 4×4 em variáveis logarítmicas, passo limitado a ±3, amortecimento em 6 níveis (1 a 0,03) e tolerância relativa 1e-11 (aceita se < 1e-8)."),
(S04, "Variáveis primárias logarítmicas", "y = ln([Zn2+], [Cl-], [NH3], [H+]) em mol/L. Mantêm as concentrações positivas e linearizam os expoentes das espécies."),
(S04, "Resíduo da especiação", "Diferença relativa entre os totais calculados (Zn, Cl, NH3~, NH4~) e os alvos. Células que não convergem tentam partida a frio; se falham, o passo é rejeitado."),
(S04, "Partida a quente e a frio", "A quente: usa y do passo anterior. A frio: estimativa a partir de NH3~ e NH4~ (NH3 livre ≈ NH3~ − 3,5·Zn~; H+ pela relação com NH4~)."),
(S04, "Jacobiana analítica", "Matriz 4×4 de ∂resíduos/∂y construída por produtos pesos·concentração·expoentes. Substitui a Jacobiana numérica do final5."),
(S04, "Sensibilidades S = dc/dQ", "Derivadas das 16 concentrações em relação às 3 quasipartículas, obtidas por função implícita (J⁻¹ vezes dT/dQ). Alimentam a matriz de difusão efetiva e o Newton da precipitação."),
(S04, "Eletroneutralidade", "Cada espécie tem carga igual a 2·nZn + nNH4~ − nCl. Por isso Cl total = 2·Zn~ + NH4~ vale identicamente."),
# ---------------- 05 ----------------
(S05, "Potencial de equilíbrio do Zn", "U_Zn = −0,762 V + (RT/2F)·ln[Zn2+]livre. Como o Zn2+ livre é minúsculo em amônia (~1e-8 M), U_Zn fica perto de −1,0 V."),
(S05, "Potencial de equilíbrio do GDE", "U_G = 1,229 V + (RT/4F)·ln(pO2) + (RT/F)·ln[H+]. Cai ≈ 59 mV por unidade de pH: pH alto no GDE reduz U_G."),
(S05, "Tensão de circuito aberto (OCV)", "Com corrente zero, V = U_G − U_Zn, médias ponderadas pela área ativa. Fica em ~1,65 a 1,7 V para os eletrólitos A a D (sem atividades iônicas)."),
(S05, "Butler-Volmer", "r = i0·A·[exp(α·η·F/RT) − exp(−α·η·F/RT)] (A/m³). Forma simétrica, com o mesmo α nos dois sentidos; aplicada ao Zn e ao GDE. No final5 a função bv() existia mas não era usada."),
(S05, "Sobretensão η", "η = φ_s − φ_e − U, onde φ_s é o potencial do eletrodo (um escalar resolvido por Newton para cada eletrodo) e φ_e o potencial local do eletrólito."),
(S05, "i0 do Zn", "2F·8e-6 = 1,54 A/m² de área ativa, de k0 do SI do artigo; α = 0,5."),
(S05, "i0 da ORR (calibrado)", "1,86e-6 A/m² de área ativa, α = 0,5. Calibrado para V ≈ 0,90 V a 1 mA/cm² no eletrólito C (leitura da Fig. 8f), não medido."),
(S05, "i0 da OER (calibrado)", "4,49e-3 A/m² de área ativa, α = 0,5. Calibrado para V ≈ 2,25 V a 2,5 mA/cm² na Célula A (leitura da Fig. 5), não medido."),
(S05, "Distribuição da corrente", "Ponto fixo em 3 iterações: pesos pela área → corrente do eletrólito → φ_e → φ_s por Newton → novas taxas. O Zn reage onde há área ativa e o GDE em toda a espessura."),
(S05, "Tensão da célula", "V = φ_s(GDE) − φ_s(Zn), incluindo equilíbrio, sobretensões e queda ôhmica/difusiva do eletrólito. Em carga a sobretensão da OER usa i0 próprio."),
(S05, "Fonte de Zn", "Dissolução: Zn~ aumenta r/(2F) por volume (r em A/m³). Na carga, r < 0 deposita Zn e eps_zn aumenta."),
(S05, "Fonte do GDE (ORR e OER)", "Descarga: NH4~ → NH3~, um mol por mol de elétrons (Tabela 5). Carga: sentido inverso. Em ambas, a corrente total do GDE é igual e oposta à do Zn."),
(S05, "ORR e OER", "Redução e evolução de O2 no GDE. O O2 não tem transporte: a pressão é constante e o fornecimento é igual à demanda."),
(S05, "Potencial local do eletrólito φ_e", "Integrado de x = 0 (φ_e = 0) pelos gradientes de potencial de cada face, incluindo o termo de difusão. A queda ôhmica total é impressa implicitamente em V."),
# ---------------- 06 ----------------
(S06, "Nernst-Planck por espécie", "Fluxo de cada espécie = difusão + migração, com o campo elétrico obtido por corrente imposta: Σ z·N = j/F. Os fluxos das quasipartículas são combinações Σ τ·N das espécies."),
(S06, "Coeficiente de difusão D", "Valor por espécie da tabela SPECIES (por exemplo Cl- 2,03e-9 e NH4+ 1,95e-9 m²/s). Dos complexos amínicos e ternários, os valores são estimados."),
(S06, "Difusão efetiva", "D vezes eps_e elevado a 1,5, média entre volumes vizinhos na face."),
(S06, "Migração e número de transferência", "Fração da corrente carregada por cada espécie, t = z²·D·c / Σ z²·D·c. A migração leva Zn²+ para o GDE e Cl- para o Zn durante a descarga."),
(S06, "Potencial de difusão", "Parcela do campo elétrico que compensa gradientes de concentração de íons com D diferentes. Está embutida no cálculo de φ_e e na matriz de difusão."),
(S06, "Condutividade do eletrólito κ", "κ = F²/(RT)·Σ z²·D·c (Nernst-Einstein), da ordem de 30 S/m para os eletrólitos A a D."),
(S06, "Matriz de difusão efetiva Dm (3×3)", "Em cada face, Dm = Σ τ·D·S − M ⊗ B, que acopla as 3 quasipartículas e inclui o potencial de difusão. É congelada no início do passo."),
(S06, "Corrente do eletrólito j_e(x)", "Soma acumulada das fontes de corrente ao longo de x. É 0 nas pontas, vale a corrente aplicada no separador e define a migração."),
(S06, "Condição de fluxo nulo", "Fluxo zero em x = 0 e x = L para Zn~, NH3~ e NH4~, incluindo a migração."),
(S06, "Sem convecção, Darcy ou Leverett", "O v2 não tem convecção, pressão capilar nem saturação. Os parâmetros Leverett e Darcy do final5 não foram portados."),
# ---------------- 07 ----------------
(S07, "Precipitação em equilíbrio local", "Em cada volume, após o transporte, cada sólido supersaturado precipita até SI = 0; sólidos presentes e subsaturados dissolvem até SI = 0 ou até acabarem. Não há cinética kr (que existia no final5)."),
(S07, "Índice de saturação SI", "SI = log10(Q) − log10(Ksp), com Q do produto iônico de espécies livres. SI > 0 precipita, SI < 0 dissolve se houver sólido."),
(S07, "Ksp e deslocamento dos hidróxidos", "Tabela 2 do artigo: ZnO −16,7; ZnCl2·4Zn(OH)2·H2O −14,2; Zn(OH)2 −17; ZnCl2·2NH3 −6,42. No código, log Ksp dos hidróxidos é somado a KSP_SHIFT_OH = 1,0 (calibração); ZnCl2·2NH3 não é alterado."),
(S07, "Produto iônico por Zn (simonkolleite)", "Q = [Zn2+]·[OH-]^1,6·[Cl-]^0,4 (5 Zn, 8 OH e 2 Cl por fórmula, normalizado por Zn). No final5 a forma completa nunca era atingida."),
(S07, "ZnO e Zn(OH)2", "Q = [Zn2+]·[OH-]². Zn(OH)2 é menos solúvel (Ksp menor) e precipita antes; ZnO só aparece se Zn(OH)2 estiver desabilitado."),
(S07, "ZnCl2·2NH3", "Q = [Zn2+]·[Cl-]²·[NH3]². Sólido cristalino com risco de passivação; consome cloreto do eletrólito."),
(S07, "ZnCl2·4Zn(OH)2·H2O (simonkolleite)", "Hidroxicloreto amorfo. Precipita em Cl alto e pH baixo, e consome cloreto e hidroxila."),
(S07, "Reações dos sólidos em Q", "Por mol de sólido: ZnO e Zn(OH)2: ΔZn~ = −1, ΔNH3~ = −2, ΔNH4~ = +2. ZnCl2·2NH3: −1, −2, 0. Simonkolleite: −5, −8, +8 (Tabela 5)."),
(S07, "Solução por Newton com intervalo", "Para cada sólido e conjunto de volumes: Newton em ξ (extensão da reação) com bissecção de reserva, até 14 iterações e |SI| < 2e-7, em até 3 varreduras sobre os sólidos."),
(S07, "Dissolução na carga", "Quando o pH ou o Zn dissolvido cai, SI < 0 e o sólido volta à solução até esgotar. É o que faz a carga recuperar Zn."),
(S07, "Sólido integrado (mol/m²) e campo (mol/m³)", "O campo é mol por m³ de volume total, por volume finito. O integral sobre dx dá mol/m² de célula, que o relatório imprime."),
(S07, "Volumes molares", "ZnO 14,5; Zn(OH)2 32,6; simonkolleite 172 e ZnCl2·2NH3 85 cm³/mol (este último é estimado). Entram no cálculo da casca de passivação."),
(S07, "Precipitado predominante", "Fase com maior quantidade integrada. 'Nenhum precipitado' se o total é < 1e-9 mol/m²."),
(S07, "Classificação do produto", "FAVORÁVEL para ZnO ou Zn(OH)2; NÃO FAVORÁVEL para fases com Cl ou NH3 (consomem cloreto e tamponam mal)."),
# ---------------- 08 ----------------
(S08, "Partícula esférica de Zn", "No modo paste cada partícula tem diâmetro inicial de 100 µm; o raio segue (eps_zn/eps_zn0)^(1/3) e a área específica A = (6/d)·eps0^(1/3)·eps^(2/3)."),
(S08, "Fração de Zn eps_zn", "Diminui na dissolução (volume molar do Zn 9,16 cm³/mol) e aumenta na deposição. Piso de 1e-3 no cálculo da área, para permitir nucleação."),
(S08, "Casca de produto e espessura δ", "Os sólidos do volume formam uma casca esférica sobre cada partícula; δ = raio externo − raio do Zn, calculado pelo volume sólido por partícula."),
(S08, "Passivação", "Fator de atividade f = 1 − (δ/δc)², limitado a [0, 1], com δc = 5 µm. Multiplica a área ativa na descarga. Quando toda a área some, a tensão colapsa e a descarga termina por 'passivação'. Na carga f = 1."),
(S08, "Área ativa do Zn em foil", "rugosidade / dx do primeiro volume do separador (rugosidade = 1). Não há aumento de área por deposição porosa; o salto de tensão entre o 1º e o 2º ciclo do artigo não é reproduzido."),
(S08, "Área específica do GDE", "4,5e3 m²/m³, uniforme na espessura do GDE."),
(S08, "GDE — eletrodo de difusão gasosa", "Região com eps_e = 0,25 e área ativa constante; o único acoplamento químico é a fonte NH4~ ⇄ NH3~. Sem inundação, sem balanço de água e sem degradação do catalisador."),
(S08, "Frente ácida na carga", "Quando o NH3 disponível no GDE acaba, o H+ da OER não é tamponado e o pH do volume cai para ~1 em segundos a minutos e a frente avança depois pelo separador. É o que limita a carga em separadores espessos."),
# ---------------- 09 ----------------
(S09, "Protocolo e segmento", "Lista de segmentos (descarga, carga ou repouso), cada um com corrente, tempo máximo, tensão de corte e, opcionalmente, capacidade máxima. A interface monta descarga + carga repetidas por 'Ciclos'."),
(S09, "Ciclos", "Repetições de descarga + carga com o mesmo estado interno. Com corrente de carga 0 só há uma descarga."),
(S09, "Tensão de corte", "Descarga termina quando V < V_corte; carga quando V > V_corte. Corte 0 V em descarga, na prática, desativa o critério (V quase nunca fica negativo)."),
(S09, "Motivo do fim", "tempo, capacidade, V_corte, passivação (Zn sem área ativa), interrompido ou max_steps."),
(S09, "Tensão final", "V do último passo do ensaio, em volts."),
(S09, "Capacidade", "Integral da corrente no tempo. O resumo soma as descargas, em Ah/m²; os segmentos mostram mAh/cm²."),
(S09, "Energia", "Integral de I·V·dt, em Wh/m². A carga entra com sinal negativo."),
(S09, "Densidade de energia (Wh/L)", "Energia da descarga dividida pela espessura total da célula (Zn + separador + GDE) e por 1000. Não inclui eletrólito consumido nem embalagem."),
(S09, "pH final do GDE", "pH médio dos volumes do GDE no fim do ensaio, a partir de [H+]."),
(S09, "pH mínimo e máximo no ensaio", "Extremos do pH médio do GDE ao longo de todos os passos."),
(S09, "Variação de pH e estabilidade", "Variação = máximo − mínimo. < 0,5: ESTÁVEL; < 1,0: MODERADAMENTE ESTÁVEL; ≥ 1,0: INSTÁVEL (mesmos limiares do final5)."),
(S09, "Perfis espaciais", "No TXT: posição, região, pH, Zn~, NH3~, NH4~ e eps_zn de cada volume no estado final."),
(S09, "Gráfico de curvas", "V (azul) e pH do GDE (laranja) contra o tempo, redesenhados a cada ~0,25 s de relógio durante a execução."),
# ---------------- 10 ----------------
(S10, "Varredura de espessura de Zn", "Descarga completa para cada espessura da lista, com o eletrólito, a corrente e o corte já preenchidos. Sempre em pasta (paste), com 4 volumes por mm (mínimo 16)."),
(S10, "Energia × espessura (Fig. 9 do artigo)", "Zn fino passiva cedo; Zn espesso sofre gradiente de pH. A tendência crescente em espessuras baixas é reproduzida; a posição do pico não foi confirmada além de 5 mm."),
(S10, "Custo da varredura", "Cada ponto leva de segundos a alguns minutos no PC e cresce com a espessura. PARAR encerra a varredura no ponto atual."),
# ---------------- 11 ----------------
(S11, "Concentrações no lugar de atividades", "Nenhum coeficiente de atividade, apesar de força iônica de 3 a 6 M. Afeta Ksp, β, Nernst e o OCV."),
(S11, "Calibrações (não são valores do artigo)", "KSP_SHIFT_OH = 1,0; i0 da ORR; i0 da OER. As duas últimas foram ajustadas a tensões lidas de figuras (Fig. 5 e Fig. 8), sem dados numéricos publicados."),
(S11, "Sem inundação, água, dendritos ou rugosidade crescente", "eps_e é estático; o inchamento por sólidos não empurra eletrólito para o GDE; não há balanço de H2O, nem morfologia dendrítica, nem aumento de área da folha de Zn."),
(S11, "Desvios conhecidos em relação ao artigo", "Na Célula A a tensão de descarga cai mais (0,61 → 0,36 V em 4 h) e a de carga sobe mais (até ~2,8 V) que na Fig. 5; o eletrólito A resulta otimista em energia (a instabilidade de pH da Fig. 8d não é reproduzida)."),
(S11, "Escopo do modelo", "Simula especiação, pH, transporte e precipitação sob corrente imposta. Não prevê vida útil, falha nem probabilidade; ciclos repetem o protocolo sobre o mesmo estado."),
# ---------------- 12 ----------------
(S12, "REF-01 — Artigo principal", "Clark, S.; Latz, A.; Horstmann, B. Rational Development of Neutral Aqueous Electrolytes for Zinc–Air Batteries. ChemSusChem, 2017, 10, 4735–4747. DOI: 10.1002/cssc.201701468. Base do modelo de eletrólito ZnCl2–NH4Cl, quasipartículas, precipitação e resultados de simulação."),
(S12, "REF-02 — Supporting Information", "Clark, S.; Latz, A.; Horstmann, B. Electronic Supporting Information: Rational Development of Neutral Aqueous Electrolytes for Zinc–Air Batteries. Equações de transporte, cinética, parâmetros físicos e implementação numérica."),
]



# ============================================================
# CLASSIFICAÇÃO — mesmos limiares do final5 (span de pH do GDE; precipitado predominante)
# ============================================================
SOLID_LABELS = ["ZnO", "ZnCl2·4Zn(OH)2·H2O", "Zn(OH)2", "ZnCl2·2NH3"]

def classify(model):
    ph0 = float(model.pH()[model.ig].mean())
    dx = np.asarray(model.dx, dtype=float)
    zn_solid = [float(np.dot(model.s[k, :], dx)) for k in range(4)]
    if max(zn_solid, default=0.0) <= 1e-9:
        main = "Nenhum precipitado"
    else:
        main = SOLID_LABELS[max(range(4), key=lambda k: zn_solid[k])]
    ph_span = model.pH_max_run - model.pH_min_run
    if ph_span < 0.5:
        phclass = "ESTÁVEL"
    elif ph_span < 1.0:
        phclass = "MODERADAMENTE ESTÁVEL"
    else:
        phclass = "INSTÁVEL"
    if main == "Nenhum precipitado":
        prod = "SEM PRECIPITADO DETECTADO"
    elif main == "ZnO":
        prod = "FAVORÁVEL (ZnO)"
    elif main == "Zn(OH)2":
        prod = "FAVORÁVEL (Zn(OH)2)"
    else:
        prod = "NÃO FAVORÁVEL (precipitado contendo Cl/NH3)"
    return phclass, prod, "CALCULADA", main, ph0, ph_span, zn_solid

S13 = "13 — VALIDAÇÃO EXPERIMENTAL (2019)"
GLOSSARY_TERMS += [
(S13, "ETAPA 2 — papel da validação", "O simulador (Etapa 1) reproduz o modelo de 2017. A Etapa 2 confronta esse mesmo modelo, sem alterar parâmetros, com os experimentos de Clark, Mainar et al. (2019). Mede onde o modelo acerta e onde diverge: os desvios são resultado, não erro de integração."),
(S13, "Artigo de 2019 (Clark, Mainar et al.)", "arXiv:1902.04400. Combina simulação e experimentos (ciclagem de célula completa, pH operando, XRD, SEM e EDS) em quatro eletrólitos de ZnCl₂–NH₄Cl. Confirma as duas previsões de 2017: pH que fica ácido na carga e produto dominante diferente de ZnO."),
(S13, "Eletrólitos E4, E6, E7 e E8", "Tabela 2 do artigo. E4: 0,51 M ZnCl₂ e 2,34 M NH₄Cl, pH 4. E6: mesma composição, pH 6 (= eletrólito A). E7: 0,26 M e 5,00 M, pH 7 (= B). E8: 0,50 M e 1,60 M, pH 8 (= C). O E4 não existia nas predefinições do simulador."),
(S13, "Células C-EI e C-OpH", "C-EI (ciclagem e XRD): eletrodos a 0,9 cm, 1,1 mL de eletrólito. C-OpH (pH operando): 2,8 cm, 4,4 mL, com microeletrodos de pH. Ambas com área ativa de 1,327 cm² e folha de Zn de 250 µm. No modelo: Zn em folha ('foil') e separador de 9 mm e 28 mm, que representam o caminho iônico entre os eletrodos."),
(S13, "Condutividade iônica (Tabela 2)", "Medida no eletrólito: 215, 206, 382 e 209 mS/cm (E4, E6, E7, E8). O modelo calcula κ a partir das concentrações das espécies e de suas difusividades (Nernst-Einstein). O teste 1 compara os valores e checa se o E7 é o mais condutivo, como medido."),
(S13, "Tabela 3 — ciclo 1", "ORR e OER (V) e sobrepotencial ΔV = OER − ORR no primeiro ciclo a 1 mA/cm², 2 h + 2 h, célula C-EI. O artigo não diz se os valores são médias, início ou fim do semiciclo."),
(S13, "Critério de tensão", "O valor do artigo deve cair dentro da faixa de tensão do modelo no semiciclo, com folga de ±0,10 V (premissa adotada porque a Tabela 3 não especifica). O critério é permissivo; por isso o placar também mostra o sobrepotencial médio, mais exigente."),
(S13, "Sobrepotencial médio ΔV", "Diferença entre a tensão média de carga e a de descarga no modelo, comparada com o ΔV da Tabela 3. A ordem entre os eletrólitos também é conferida."),
(S13, "pH operando (Figs. 7 e 8)", "pH medido por microeletrodos perto do eletrodo de ar e do Zn durante um ciclo a 2 mA/cm², 29 h de descarga e 29 h de carga, célula C-OpH. A medição é pontual, e o modelo é 1D: o artigo afirma que o modelo tende a superestimar as variações."),
(S13, "Direções de pH (Fig. 5)", "No eletrólito E8 a 1 mA/cm²: o pH do ar alcaliniza na descarga e acidifica na carga; o do Zn acidifica na descarga e alcaliniza na carga. O teste 2 confere o sentido de cada variação."),
(S13, "Platô de pH do E4", "O E4 parte de pH 4 e, no experimento, sobe rápido no início da descarga e se aproxima de um valor estável perto de pH 6 (primeiro no ar, depois no Zn). Critério: pH ao fim da descarga entre 5,5 e 6,5 (premissa)."),
(S13, "Queda ácida na carga", "No E4 o pH cai a valores fortemente ácidos perto do fim da carga, começando no ar e seguindo no Zn. Segundo o artigo, a causa é o transporte lento e o esgotamento local de NH₃. Critério: queda de pelo menos 1,0 unidade de pH (limiar de INSTÁVEL de classify())."),
(S13, "DoD de 40% (58 mAh/cm²)", "Protocolo do XRD: descarga a 1 mA/cm² até 40% de profundidade de descarga. Equivale a 58 mAh/cm², igual a 2 mA/cm² × 29 h do ensaio operando. Premissa: célula C-EI, que o texto não especifica para o XRD."),
(S13, "Proxy de XRD", "O XRD não é quantitativo. O teste 4 usa a fração molar de Zn presa em cada sólido para decidir qual fase domina e compara com a fase dominante do artigo."),
(S13, "Fase em camadas intercalada (E8)", "No E8 o artigo observa uma fase em camadas com espaçamento c = 10,68 Å, compatível com hidróxido-cloreto de zinco com NH₃ intercalada, não identificada nos bancos de dados. Ela não existe no conjunto de 4 sólidos do modelo: desvio conhecido (status CONHECIDO)."),
(S13, "Ausência de ZnO", "O artigo não identifica ZnO em nenhuma amostra descarregada. O teste exige que a fração de Zn em ZnO seja menor que 5%."),
(S13, "Autofreio (Fig. 12)", "Ao depositar Zn a partir dos complexos Zn(NH₃)ₓ, o NH₃ é liberado, mantém o pH perto da neutralidade e a solubilidade do Zn baixa, o que retarda a redissolução dos precipitados na carga."),
(S13, "Retenção e perda líquida de Zn", "No experimento, após a carga resta precipitado e o Zn dissolvido fica abaixo do inicial. O teste 5 mede a fração de precipitado que sobra e compara o Zn dissolvido final com o inicial. As condições da Fig. 12 não constam do texto: usa-se E6, C-EI, 1 mA/cm², 10 h + 10 h."),
(S13, "Balanço de Zn (Faraday)", "Zn dissolvido + Zn em sólidos deve igualar o inicial + capacidade líquida/(2F). Erro acima de 2% do Zn ciclado indica perda de conservação numérica."),
(S13, "Esgotamento de Zn no eletrodo de Zn", "Se o Zn~ do lado do Zn chega ao piso numérico (ZT_FLOOR), o modelo segue impondo a corrente de deposição, não há outra reação e o balanço de Zn deixa de fechar. Observado na carga de 29 h do E4. Marca o limite de validade da carga longa; o piso numérico é a causa provável."),
(S13, "Status do placar", "OK: reproduz o artigo dentro do critério. DIVERGE: não reproduz. CONHECIDO: desvio já identificado. INDICATIVO: modo reduzido, sem valor de comparação. PENDENTE: não executado ou condição pouco informativa. INFO: observação sem nota."),
(S13, "Modo completo e reduzido", "Completo usa as condições do artigo (29 h + 29 h; 40% de DoD; 10 h + 10 h). Reduzido encurta esses três testes para rodar rápido; o resultado vira INDICATIVO. Os testes 1 e 2 são sempre completos."),
(S13, "Por que não recalibrar", "Ajustar i0 e demais parâmetros nos dados da Tabela 3 e depois validar nos mesmos dados seria circular. A Etapa 2 mostra o modelo como está; uma recalibração, se feita, fica como etapa própria."),
(S13, "Fora do alcance", "Evaporação do eletrólito (até 14% no E8, Tabela 3), inundação do eletrodo de ar, dissolução do MnO₂ em meio ácido, reação de evolução de H₂ e degradação ao longo de 150 ciclos. O O₂ dissolvido medido (6,2 a 7,1 mg/L) também não entra: o modelo usa pO₂ constante."),
(S12, "REF-03 — Artigo experimental (2019)", "Clark, S.; Mainar, A. R.; Iruin, E.; Colmenares, L. C.; Blázquez, J. A.; Tolchard, J. R.; Latz, A.; Horstmann, B. Towards Rechargeable Zinc-Air Batteries with Aqueous Chloride Electrolytes. arXiv:1902.04400, fev. 2019. Base da Etapa 2. Material suplementar não consultado."),
]

# ============================================================
# PRESETS (valores do artigo; protocolo = o usado nas figuras)
# ============================================================
PRESETS = {
    "Eletrólito A — célula padrão": {"ZnCl₂ (mol/L)": 0.51, "NH₄Cl (mol/L)": 2.34, "pH inicial": 6.0, "Modelo de Zn": "paste", "Sólidos": "todos",
                                      "Zn (mm)": 5.0, "Separador (mm)": 0.5, "GDE (mm)": 0.5, "Ciclos": 1,
                                      "I descarga (mA/cm²)": 1.0, "t máx descarga (h)": 2000, "V corte descarga (V)": 0.6,
                                      "I carga (mA/cm²)": 0.0, "t máx carga (h)": 8, "V corte carga (V)": 2.6},
    "Eletrólito B — célula padrão": {"ZnCl₂ (mol/L)": 0.26, "NH₄Cl (mol/L)": 5.0, "pH inicial": 7.0, "Modelo de Zn": "paste", "Sólidos": "todos",
                                      "Zn (mm)": 5.0, "Separador (mm)": 0.5, "GDE (mm)": 0.5, "Ciclos": 1,
                                      "I descarga (mA/cm²)": 1.0, "t máx descarga (h)": 2000, "V corte descarga (V)": 0.6,
                                      "I carga (mA/cm²)": 0.0, "t máx carga (h)": 8, "V corte carga (V)": 2.6},
    "Eletrólito C — célula padrão": {"ZnCl₂ (mol/L)": 0.5, "NH₄Cl (mol/L)": 1.6, "pH inicial": 8.0, "Modelo de Zn": "paste", "Sólidos": "todos",
                                      "Zn (mm)": 5.0, "Separador (mm)": 0.5, "GDE (mm)": 0.5, "Ciclos": 1,
                                      "I descarga (mA/cm²)": 1.0, "t máx descarga (h)": 2000, "V corte descarga (V)": 0.6,
                                      "I carga (mA/cm²)": 0.0, "t máx carga (h)": 8, "V corte carga (V)": 2.6},
    "Eletrólito D — só ZnO (ideal)": {"ZnCl₂ (mol/L)": 0.5, "NH₄Cl (mol/L)": 2.0, "pH inicial": 7.0, "Modelo de Zn": "paste", "Sólidos": "só ZnO",
                                      "Zn (mm)": 5.0, "Separador (mm)": 0.5, "GDE (mm)": 0.5, "Ciclos": 1,
                                      "I descarga (mA/cm²)": 1.0, "t máx descarga (h)": 2000, "V corte descarga (V)": 0.6,
                                      "I carga (mA/cm²)": 0.0, "t máx carga (h)": 8, "V corte carga (V)": 2.6},
    "Célula A do artigo — 3 ciclos": {"ZnCl₂ (mol/L)": 0.51, "NH₄Cl (mol/L)": 2.34, "pH inicial": 6.0, "Modelo de Zn": "foil", "Sólidos": "todos",
                                      "Zn (mm)": 0.5, "Separador (mm)": 30.0, "GDE (mm)": 0.5, "Ciclos": 3,
                                      "I descarga (mA/cm²)": 5.0, "t máx descarga (h)": 4, "V corte descarga (V)": 0.0,
                                      "I carga (mA/cm²)": 2.5, "t máx carga (h)": 8, "V corte carga (V)": 4.0},
}


# ============================================================
# ETAPA 2 — VALIDAÇÃO EXPERIMENTAL (2019)
# Clark, Mainar et al., arXiv:1902.04400. Mesmo modelo v2; nenhuma calibração nova.
# Parte numérica pura (sem Tkinter): condições, critérios, execução, placar e relatório.
# ============================================================
VAL_ELETROLITOS = {"E4": (0.51, 2.34, 4.0), "E6": (0.51, 2.34, 6.0),      # Tabela 2: ZnCl2 (M), NH4Cl (M), pH
                   "E7": (0.26, 5.00, 7.0), "E8": (0.50, 1.60, 8.0)}
VAL_LABEL = {"E4": "E4", "E6": "E6 = A", "E7": "E7 = B", "E8": "E8 = C"}
VAL_COR = {"E4": "#c62828", "E6": "#1976d2", "E7": "#2e7d32", "E8": "#ef6c00"}
VAL_GEOM = {"C-EI": 9.0, "C-OpH": 28.0}                       # distância entre eletrodos (mm): 0,9 cm e 2,8 cm; folha de Zn
VAL_IC = {"E4": 215.0, "E6": 206.0, "E7": 382.0, "E8": 209.0}  # Tabela 2: condutividade iônica medida (mS/cm)
VAL_TAB3 = {"E4": (0.887, 2.079), "E6": (0.929, 2.067),       # Tabela 3, ciclo 1: (ORR, OER) em V, 1 mA/cm2, 2 h + 2 h
            "E7": (0.860, 1.977), "E8": (0.888, 1.922)}
VAL_XRD = {"E4": "ZnCl2_4ZnOH2_H2O", "E6": "ZnCl2_4ZnOH2_H2O", "E7": "ZnCl2_2NH3", "E8": None}   # Fig. 9: fase dominante
VAL_XRD_TXT = {"E4": "simonkolleíta", "E6": "simonkolleíta", "E7": "ZnCl₂·2NH₃",
               "E8": "camada intercalada com NH₃ (fora do modelo)"}
VAL_FASE_NOME = {"ZnO": "ZnO", "ZnCl2_4ZnOH2_H2O": "simonkolleíta", "ZnOH2": "Zn(OH)₂", "ZnCl2_2NH3": "ZnCl₂·2NH₃"}
VAL_FASE_COR = ["#6a1b9a", "#c62828", "#1976d2", "#212121"]   # mesma ordem de SOLIDS
VAL_ZN_POR_SOLIDO = np.array([1.0, 5.0, 1.0, 1.0])            # mol de Zn por mol de sólido (ordem de SOLIDS)
VAL_ASSUME = {
    "tol_kappa": 0.25,             # tolerância relativa da condutividade vs Tabela 2
    "tol_V": 0.10,                 # V: tolerância sobre a faixa de tensão do modelo vs Tabela 3
    "plato_pH": (5.5, 6.5),        # pH do E4 ao fim da descarga ("perto de pH 6")
    "limiar_pH": 1.0,              # mesmo limiar de INSTÁVEL de classify()
    "zno_max": 0.05,               # fração de Zn em ZnO abaixo da qual "não há ZnO"
    "dod40_mAh_cm2": 58.0,         # 40% de DoD = 2 mA/cm2 x 29 h do ensaio operando
    "faraday_tol": 0.02,           # erro do balanço de Zn, em fração do Zn ciclado
    "zn_esgotado": 0.01,           # mol/m3: Zn~ no lado do Zn abaixo disso = esgotado
}
VAL_TESTS = [("t1", "1 · CONDUTIVIDADE IÔNICA (Tabela 2)"),
             ("t2", "2 · TENSÕES E DIREÇÕES DE pH (Tabela 3, Fig. 5)"),
             ("t3", "3 · pH OPERANDO (Figs. 7 e 8)"),
             ("t4", "4 · FASES SÓLIDAS A 40% DE DoD (Fig. 9)"),
             ("t5", "5 · AUTOFREIO E BALANÇO DE Zn (Fig. 12)")]
VAL_STATUS_COR = {"OK": "#2e7d32", "DIVERGE": "#c62828", "CONHECIDO": "#ef6c00",
                  "INDICATIVO": "#546e7a", "PENDENTE": "#9e9e9e", "INFO": "#1565c0"}


def val_hours(tid, modo):
    """Horas simuladas de cada teste (peso da barra de progresso)."""
    red = modo == "reduzido"
    return {"t1": 0.05, "t2": 16.0, "t3": 8.0*(6.0 if red else 29.0),
            "t4": 4.0*(15.0 if red else VAL_ASSUME["dod40_mAh_cm2"]), "t5": 2.0*(4.0 if red else 10.0)}[tid]


class ValCtl:
    """Controle da validação: parada, progresso por horas simuladas e menor Zn~ visto no lado do Zn."""
    def __init__(self):
        self.stop = False
        self.on_progress = None
        self.on_partial = None
        self.total_h = 1.0
        self.done_h = 0.0
        self.zmin = float("inf")

    def tick(self, unit_h, frac, label, V=None):
        if self.on_progress:
            self.on_progress(min(1.0, (self.done_h + unit_h*frac)/max(self.total_h, 1e-9)), label, V)

    def partial(self, tid, r):
        if self.on_partial:
            self.on_partial(tid, r)


def val_new_state():
    return {"res": {}, "checks": [], "modos": {}, "running": False, "tests": [], "elapsed": 0.0, "report": ""}


def val_model(el, geom):
    zn, nh, ph = VAL_ELETROLITOS[el]
    return Model({"zncl2_M": zn, "nh4cl_M": nh, "initial_pH": ph, "zn_model": "foil", "zn_mm": 0.5,
                  "separator_mm": VAL_GEOM[geom], "gde_mm": 0.5})


def val_seg(mode, I, h, vcut=None):
    return {"mode": mode, "I_mA_cm2": I, "t_max_h": h,
            "V_cut": vcut if vcut is not None else (0.0 if mode == "discharge" else 4.0)}   # sem corte prematuro


def val_zn_diss(m):        # mol/m2 de Zn dissolvido
    return float(np.dot(m.dx, m.eps_e*m.Q[:, 0]))


def val_zn_sol(m):         # mol/m2 de Zn em cada sólido (ordem de SOLIDS)
    return np.array([np.dot(m.dx, m.s[i]) for i in range(4)])*VAL_ZN_POR_SOLIDO


def val_kappa(m):          # mS/cm, eletrólito do seio da solução (Nernst-Einstein com as espécies do modelo)
    return float(F*(m.c[m.ig[0]] @ Z2D)/RTF*10.0)


def val_pH(m):
    ph = m.pH()
    return float(ph[m.ig].mean()), float(ph[0])


def val_stats(h):
    v = np.array([x["V"] for x in h])
    return {"ini": float(v[len(v)//10]), "med": float(v.mean()), "fim": float(v[-1]),
            "vmin": float(v.min()), "vmax": float(v.max())}


def val_ds(rows, maxn=240):
    if len(rows) <= maxn:
        return list(rows)
    k = int(math.ceil(len(rows)/maxn))
    out = rows[::k]
    return list(out) + ([rows[-1]] if out[-1] is not rows[-1] else [])


def val_run(m, segs, ctl, unit_h, label, rec=60.0):
    def cb(mm, k, n, frac, V, cpu):
        ctl.tick(unit_h, (k + frac)/n, label, V)
        z = float(mm.Q[0, 0])
        if z < ctl.zmin:
            ctl.zmin = z
        if ctl.stop:
            mm.stop_flag = True
    out = m.run_protocol(segs, rec_every_s=rec, callback=cb)
    ctl.done_h += unit_h
    return out


def val_chk(out, tid, txt, st, det="", red=False):
    if red and st in ("OK", "DIVERGE"):
        det = f"[modo reduzido: resultado {st}] " + det
        st = "INDICATIVO"
    out.append({"t": tid, "txt": txt, "st": st, "det": det})


def val_t1(ctl, modo, out, res):
    r = {el: (val_kappa(val_model(el, "C-EI")), VAL_IC[el]) for el in VAL_ELETROLITOS}
    res["t1"] = r
    ctl.done_h += val_hours("t1", modo)
    tol = VAL_ASSUME["tol_kappa"]
    for el, (km, ka) in r.items():
        val_chk(out, "t1", f"{VAL_LABEL[el]}: κ modelo {km:.0f} contra {ka:.0f} mS/cm",
                "OK" if abs(km/ka - 1) <= tol else "DIVERGE", f"razão {km/ka:.2f} (tolerância ±{100*tol:.0f}%)")
    top_m, top_a = max(r, key=lambda e: r[e][0]), max(r, key=lambda e: r[e][1])
    val_chk(out, "t1", "E7 é o eletrólito mais condutivo", "OK" if top_m == top_a == "E7" else "DIVERGE",
            f"modelo: {top_m}; artigo: {top_a}")


def val_t2(ctl, modo, out, res):
    tol = VAL_ASSUME["tol_V"]; r = {}
    for el in VAL_ELETROLITOS:
        if ctl.stop:
            break
        m = val_model(el, "C-EI"); g0, z0 = val_pH(m)
        val_run(m, [val_seg("discharge", 1.0, 2.0)], ctl, 2.0, f"T2 · {el} descarga", rec=30.0)
        if ctl.stop:
            break
        nd = len(m.history); hd = list(m.history); gd, zd = val_pH(m)
        val_run(m, [val_seg("charge", 1.0, 2.0)], ctl, 2.0, f"T2 · {el} carga", rec=30.0)
        if ctl.stop:
            break
        gc, zc = val_pH(m)
        r[el] = {"g0": g0, "z0": z0, "gd": gd, "zd": zd, "gc": gc, "zc": zc,
                 "D": val_stats(hd), "C": val_stats(m.history[nd:]),
                 "hist": val_ds([(h["t_h"], h["V"], h["pH_gde"], h["pH_zn"]) for h in m.history])}
        res["t2"] = r; ctl.partial("t2", r)
    res["t2"] = r
    for el, x in r.items():
        orr, oer = VAL_TAB3[el]; D, C = x["D"], x["C"]
        val_chk(out, "t2", f"{VAL_LABEL[el]}: ORR {orr:.3f} V na faixa do modelo",
                "OK" if D["vmin"] - tol <= orr <= D["vmax"] + tol else "DIVERGE",
                f"modelo {D['vmin']:.3f}–{D['vmax']:.3f} V (média {D['med']:.3f}) ± {tol:.2f}")
        val_chk(out, "t2", f"{VAL_LABEL[el]}: OER {oer:.3f} V na faixa do modelo",
                "OK" if C["vmin"] - tol <= oer <= C["vmax"] + tol else "DIVERGE",
                f"modelo {C['vmin']:.3f}–{C['vmax']:.3f} V (média {C['med']:.3f}) ± {tol:.2f}")
    if len(r) == len(VAL_ELETROLITOS):
        dv = lambda e: r[e]["C"]["med"] - r[e]["D"]["med"]
        da = lambda e: VAL_TAB3[e][1] - VAL_TAB3[e][0]
        om, oa = sorted(r, key=dv), sorted(r, key=da)
        val_chk(out, "t2", "ordem de sobrepotencial entre eletrólitos", "OK" if om == oa else "DIVERGE",
                "modelo " + " < ".join(om) + "; artigo " + " < ".join(oa))
        val_chk(out, "t2", "sobrepotencial médio ΔV (modelo contra artigo)", "INFO",
                "; ".join(f"{e} {dv(e):.2f} contra {da(e):.2f} V" for e in r))
    if "E8" in r:
        x = r["E8"]
        for txt, ok, det in (
            ("E8: pH do ar sobe na descarga", x["gd"] > x["g0"], f"{x['g0']:.2f} → {x['gd']:.2f}"),
            ("E8: pH do ar cai na carga", x["gc"] < x["gd"], f"{x['gd']:.2f} → {x['gc']:.2f}"),
            ("E8: pH do Zn cai na descarga", x["zd"] < x["z0"], f"{x['z0']:.2f} → {x['zd']:.2f}"),
            ("E8: pH do Zn sobe na carga", x["zc"] > x["zd"], f"{x['zd']:.2f} → {x['zc']:.2f}")):
            val_chk(out, "t2", txt, "OK" if ok else "DIVERGE", det)
    if len(r) < len(VAL_ELETROLITOS):
        val_chk(out, "t2", "teste incompleto", "PENDENTE", f"{len(r)} de {len(VAL_ELETROLITOS)} eletrólitos concluídos")


def val_t3(ctl, modo, out, res):
    red = modo == "reduzido"; hh = 6.0 if red else 29.0
    lim = VAL_ASSUME["limiar_pH"]; lo, hi = VAL_ASSUME["plato_pH"]; r = {}
    for el in VAL_ELETROLITOS:
        if ctl.stop:
            break
        m = val_model(el, "C-OpH"); g0, z0 = val_pH(m)
        tot0 = val_zn_diss(m) + float(val_zn_sol(m).sum())
        val_run(m, [val_seg("discharge", 2.0, hh)], ctl, hh, f"T3 · {el} descarga", rec=120.0)
        if ctl.stop:
            break
        nd = len(m.history); gd, zd = val_pH(m)
        ctl.zmin = float("inf")
        val_run(m, [val_seg("charge", 2.0, hh)], ctl, hh, f"T3 · {el} carga", rec=120.0)
        if ctl.stop:
            break
        gc, zc = val_pH(m); hc = m.history[nd:]
        allg = [g0] + [h["pH_gde"] for h in m.history]
        g2 = next((h["pH_gde"] for h in m.history if h["t_h"] >= 2.0), gd)
        dep = 2.0*10.0*hh*3600.0/(2*F)
        esp = tot0 + m.cap*3600.0/(2*F)
        err = abs(val_zn_diss(m) + float(val_zn_sol(m).sum()) - esp)/max(dep, 1e-12)
        r[el] = {"g0": g0, "gd": gd, "gmin": min(h["pH_gde"] for h in hc), "gc": gc, "z0": z0, "zd": zd, "zc": zc,
                 "faixa": max(allg) - min(allg), "subida2h": g2 - g0, "hh": hh,
                 "queda": gd - min(h["pH_gde"] for h in hc), "bal_err": err, "zmin": min(ctl.zmin, float(m.Q[0, 0])),
                 "hist": val_ds([(h["t_h"], h["pH_gde"], h["pH_zn"]) for h in m.history])}
        res["t3"] = r; ctl.partial("t3", r)
    res["t3"] = r
    if "E4" in r:
        e4 = r["E4"]
        val_chk(out, "t3", f"E4: pH do ar sobe ≥ {lim:.1f} e termina em {lo}–{hi}",
                "OK" if (e4["gd"] - e4["g0"] >= lim and lo <= e4["gd"] <= hi) else "DIVERGE",
                f"{e4['g0']:.2f} → {e4['gd']:.2f}", red)
        val_chk(out, "t3", f"E4: pH do Zn termina em {lo}–{hi} na descarga", "OK" if lo <= e4["zd"] <= hi else "DIVERGE",
                f"{e4['z0']:.2f} → {e4['zd']:.2f}", red)
        val_chk(out, "t3", f"E4: queda ácida do ar na carga (≥ {lim:.1f})", "OK" if e4["queda"] >= lim else "DIVERGE",
                f"queda {e4['queda']:.2f}; mínimo {e4['gmin']:.2f}", red)
        for el in ("E6", "E7", "E8"):
            if el in r:
                x = r[el]
                ok = x["faixa"] < e4["faixa"] and x["queda"] < lim and x["subida2h"] < lim
                val_chk(out, "t3", f"{VAL_LABEL[el]}: pH mais estável que o E4, sem queda ácida nem salto inicial",
                        "OK" if ok else "DIVERGE",
                        f"faixa {x['faixa']:.2f} (E4 {e4['faixa']:.2f}); queda {x['queda']:.2f}; salto 2 h {x['subida2h']:+.2f}", red)
    for el, x in r.items():
        val_chk(out, "t3", f"{VAL_LABEL[el]}: balanço de Zn fecha (Faraday)", "OK" if x["bal_err"] <= VAL_ASSUME["faraday_tol"] else "DIVERGE",
                f"erro {100*x['bal_err']:.1f}% do Zn ciclado (tolerância {100*VAL_ASSUME['faraday_tol']:.0f}%)", red)
        if x["zmin"] <= VAL_ASSUME["zn_esgotado"]:
            val_chk(out, "t3", f"{VAL_LABEL[el]}: Zn~ esgota no lado do Zn durante a carga", "DIVERGE",
                    f"mínimo {x['zmin']:.3g} mol/m³: limite do modelo (sem outra reação além do Zn)", red)
    if r:
        val_chk(out, "t3", "variações de pH do modelo contra as medidas", "INFO",
                "o artigo diz que o modelo 1D tende a superestimar as variações (medição pontual, 3D)")
    if len(r) < len(VAL_ELETROLITOS):
        val_chk(out, "t3", "teste incompleto", "PENDENTE", f"{len(r)} de {len(VAL_ELETROLITOS)} eletrólitos concluídos")


def val_t4(ctl, modo, out, res):
    red = modo == "reduzido"; q = 15.0 if red else VAL_ASSUME["dod40_mAh_cm2"]; r = {}
    for el in VAL_ELETROLITOS:
        if ctl.stop:
            break
        m = val_model(el, "C-EI")
        val_run(m, [val_seg("discharge", 1.0, q, vcut=-9.0)], ctl, q, f"T4 · {el} descarga", rec=600.0)
        if ctl.stop:
            break
        zs = val_zn_sol(m); tot = float(zs.sum())
        fr = [float(v)/tot for v in zs] if tot > 1e-9 else [0.0]*4
        r[el] = {"f": fr, "dom": SOLIDS[int(np.argmax(fr))] if tot > 1e-9 else None, "q": q}
        res["t4"] = r; ctl.partial("t4", r)
    res["t4"] = r
    for el in ("E4", "E6", "E7"):
        if el in r:
            ok = r[el]["dom"] == VAL_XRD[el]
            nome = VAL_FASE_NOME.get(r[el]["dom"], "sem sólido")
            val_chk(out, "t4", f"{VAL_LABEL[el]}: fase dominante {VAL_XRD_TXT[el]}", "OK" if ok else "DIVERGE",
                    f"modelo: {nome} ({100*max(r[el]['f']):.0f}% do Zn)", red)
    if "E8" in r:
        e = r["E8"]
        val_chk(out, "t4", "E8 = C: fase em camadas intercalada com NH₃", "CONHECIDO",
                f"modelo: {VAL_FASE_NOME.get(e['dom'], 'sem sólido')} ({100*max(e['f']):.0f}% do Zn); a fase do artigo não existe nos 4 sólidos do modelo")
    if r:
        zno = max(v["f"][0] for v in r.values())
        val_chk(out, "t4", "sem ZnO identificável em nenhuma amostra", "OK" if zno < VAL_ASSUME["zno_max"] else "DIVERGE",
                f"maior fração de Zn em ZnO: {100*zno:.1f}% (limite {100*VAL_ASSUME['zno_max']:.0f}%)", red)
    if len(r) < len(VAL_ELETROLITOS):
        val_chk(out, "t4", "teste incompleto", "PENDENTE", f"{len(r)} de {len(VAL_ELETROLITOS)} eletrólitos concluídos")


def val_t5(ctl, modo, out, res):
    red = modo == "reduzido"; hh = 4 if red else 10
    m = val_model("E6", "C-EI"); g, z = val_pH(m)
    z0 = val_zn_diss(m); tot0 = z0 + float(val_zn_sol(m).sum())
    smp = [(0.0, z0, float(val_zn_sol(m).sum()), g, z)]
    ok_run = True
    for modo_seg, nome in (("discharge", "descarga"), ("charge", "carga")):
        for k in range(hh):
            val_run(m, [val_seg(modo_seg, 1.0, 1.0)], ctl, 1.0, f"T5 · E6 {nome} {k+1}/{hh} h", rec=60.0)
            if ctl.stop:
                ok_run = False
                break
            g, z = val_pH(m)
            smp.append((m.t/3600.0, val_zn_diss(m), float(val_zn_sol(m).sum()), g, z))
            ctl.partial("t5", {"samples": list(smp)})
        if modo_seg == "discharge":
            idx_d = len(smp) - 1
        if not ok_run:
            break
    r = {"samples": smp, "hh": hh}
    res["t5"] = r
    if not ok_run:
        val_chk(out, "t5", "teste incompleto", "PENDENTE", "interrompido antes do fim da carga")
        return
    dep = 10.0*hh*3600.0/(2*F)
    _, zd, sd, _, pzd = smp[idx_d]; _, zc, sc, _, pzc = smp[-1]
    err = abs(zc + sc - (tot0 + m.cap*3600.0/(2*F)))/max(dep, 1e-12)
    r.update({"z0": z0, "zd": zd, "zc": zc, "sd": sd, "sc": sc, "dep": dep, "err": err, "pzd": pzd, "pzc": pzc})
    val_chk(out, "t5", "balanço de Zn fecha (Faraday)", "OK" if err <= VAL_ASSUME["faraday_tol"] else "DIVERGE",
            f"erro {100*err:.2f}% do Zn ciclado", red)
    if sd < 1e-3:
        val_chk(out, "t5", "houve precipitação na descarga", "PENDENTE",
                f"só {sd:.4f} mol/m² precipitaram: condição pouco informativa para o autofreio", red)
        return
    val_chk(out, "t5", "precipitado não se redissolve por completo (retenção > 0)", "OK" if sc > 1e-3 else "DIVERGE",
            f"sólido: {sd:.3f} → {sc:.3f} mol/m²; retenção {100*sc/sd:.1f}% do formado", red)
    val_chk(out, "t5", "perda líquida de Zn dissolvido ao fim da carga", "OK" if zc < z0*(1 - 1e-3) else "DIVERGE",
            f"dissolvido: início {z0:.3f}, fim da descarga {zd:.3f}, fim da carga {zc:.3f} mol/m²", red)
    val_chk(out, "t5", "autofreio: pH do Zn sobe na carga", "OK" if pzc > pzd else "DIVERGE",
            f"{pzd:.2f} → {pzc:.2f}", red)


VAL_FUNCS = {"t1": val_t1, "t2": val_t2, "t3": val_t3, "t4": val_t4, "t5": val_t5}


def val_counts(checks):
    c = {k: 0 for k in VAL_STATUS_COR}
    for x in checks:
        c[x["st"]] = c.get(x["st"], 0) + 1
    return c


def val_report(st):
    L = []

    def P(x=""):
        L.append(x)

    P("VALIDAÇÃO EXPERIMENTAL (2019) — ZAB NEUTRA (v2)")
    P("Clark, Mainar et al., arXiv:1902.04400. Mesmo modelo, sem recalibração.")
    P("=" * 60)
    P("Premissas: " + "; ".join(f"{k} = {v}" for k, v in VAL_ASSUME.items()))
    res = st["res"]
    for tid, title in VAL_TESTS:
        P(); P(title); P("-" * 60)
        r = res.get(tid)
        if not r:
            P("  (não executado)"); continue
        P(f"  modo: {st['modos'].get(tid, '—')}")
        if tid == "t1":
            for el, (km, ka) in r.items():
                P(f"  {VAL_LABEL[el]:<7} κ modelo {km:7.1f}  artigo {ka:6.1f} mS/cm  razão {km/ka:.2f}")
        elif tid == "t2":
            P("  El.     ORR art |  modelo ini   méd   fim | OER art |  modelo ini   méd   fim")
            for el, x in r.items():
                D, C = x["D"], x["C"]
                P(f"  {VAL_LABEL[el]:<7} {VAL_TAB3[el][0]:6.3f} |         {D['ini']:5.3f} {D['med']:5.3f} {D['fim']:5.3f} | "
                  f"{VAL_TAB3[el][1]:6.3f} |         {C['ini']:5.3f} {C['med']:5.3f} {C['fim']:5.3f}")
            if "E8" in r:
                x = r["E8"]
                P(f"  E8 pH ar: {x['g0']:.2f} → {x['gd']:.2f} (descarga) → {x['gc']:.2f} (carga); pH Zn: {x['z0']:.2f} → {x['zd']:.2f} → {x['zc']:.2f}")
        elif tid == "t3":
            P("  El.     pH ar: 0   fim D  mín C  fim C | Zn: fim D fim C | faixa ar | erro Zn | Zn~ mín")
            for el, x in r.items():
                P(f"  {VAL_LABEL[el]:<7} {x['g0']:11.2f} {x['gd']:6.2f} {x['gmin']:6.2f} {x['gc']:6.2f} | {x['zd']:11.2f} {x['zc']:5.2f} | "
                  f"{x['faixa']:8.2f} | {100*x['bal_err']:5.1f}% | {x['zmin']:.3g}")
        elif tid == "t4":
            P("  El.     ZnO  simonk. Zn(OH)2 ZnCl2·2NH3 | dominante (modelo) | artigo")
            for el, x in r.items():
                f = x["f"]
                P(f"  {VAL_LABEL[el]:<7} {f[0]:4.2f} {f[1]:7.2f} {f[2]:7.2f} {f[3]:10.2f} | {VAL_FASE_NOME.get(x['dom'], '—'):<17} | {VAL_XRD_TXT[el]}")
        elif tid == "t5" and "z0" in r:
            P(f"  Zn dissolvido: início {r['z0']:.3f}, fim da descarga {r['zd']:.3f}, fim da carga {r['zc']:.3f} mol/m²")
            P(f"  Zn em sólidos: fim da descarga {r['sd']:.3f}, fim da carga {r['sc']:.3f} mol/m² (Faraday da descarga: {r['dep']:.3f})")
    P(); P("PLACAR"); P("=" * 60)
    cc = val_counts(st["checks"])
    P("  " + "  ".join(f"{k} {v}" for k, v in cc.items() if v))
    for tid, title in VAL_TESTS:
        xs = [x for x in st["checks"] if x["t"] == tid]
        if not xs:
            continue
        P(); P(title)
        for x in xs:
            P(f"  [{x['st']:<10}] {x['txt']}")
            if x["det"]:
                P(f"               {x['det']}")
    P(); P("COMO LER"); P("-" * 60)
    P("  OK: o modelo reproduz o que o artigo reporta, dentro do critério declarado.")
    P("  DIVERGE: o modelo não reproduz; é resultado, não erro da integração.")
    P("  CONHECIDO: desvio já identificado (fase fora do conjunto de sólidos).")
    P("  INDICATIVO: modo reduzido; condição mais curta que a do artigo, sem valor de comparação.")
    P("  PENDENTE: não executado ou condição pouco informativa. INFO: sem nota.")
    P("  Sem recalibração: ajustar nos mesmos dados e depois validar seria circular.")
    return "\n".join(L)


class App:
    def __init__(self, root):
        self.root = root
        root.title("ZAB — Modelo v2")
        root.geometry("360x640")
        root.minsize(320, 520)

        self.running = False
        self.model = None
        self.vars = {}
        self.summary_vars = {}
        self.val_state = val_new_state()
        self.val_ctl = None

        # Paleta simples e leve para Tkinter/Android
        bg = "#f2f2f2"
        card = "#ffffff"
        root.configure(bg=bg)

        outer = tk.Frame(root, bg=bg)
        outer.pack(fill="both", expand=True)

        self.canvas = tk.Canvas(outer, bg=bg, highlightthickness=0, bd=0)
        self.canvas.pack(fill="both", expand=True)

        self.frame = tk.Frame(self.canvas, bg=bg)
        self.window_id = self.canvas.create_window((0, 0), window=self.frame, anchor="nw")
        self.frame.bind("<Configure>", self._on_frame_configure)
        self.canvas.bind("<Configure>", self._on_canvas_configure)

        # Rolagem por toque / mouse
        self._bind_touch_scroll(self.canvas)

        self._header()
        self._stages_card()
        self._input_card(card)
        self._protocol_card(card)
        self._run_card(card)
        self._summary_card(card)
        self._charts_cards(card)
        self._result_card(card)

    # ---------- layout (mesmo padrão do final5) ----------
    def _on_frame_configure(self, event=None):
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))

    def _on_canvas_configure(self, event):
        self.canvas.itemconfigure(self.window_id, width=max(300, event.width))

    def _bind_touch_scroll(self, canvas):
        """Rolagem por gesto vertical com o dedo, sem depender de scrollbar."""
        state = {"dragging": False, "last_y": 0}

        def inside_canvas(widget):
            while widget is not None:
                if widget is canvas:
                    return True
                try:
                    widget = widget.master
                except Exception:
                    return False
            return False

        def press(event):
            if not inside_canvas(event.widget):
                return
            state["dragging"] = True
            state["last_y"] = event.y_root

        def motion(event):
            if not state["dragging"]:
                return
            dy = event.y_root - state["last_y"]
            state["last_y"] = event.y_root
            if dy:
                top, bottom = canvas.yview()
                span = max(1.0, bottom - top)
                canvas.yview_moveto(top - (dy / max(1, canvas.winfo_height())) * span)

        def release(event):
            state["dragging"] = False

        canvas.bind_all("<ButtonPress-1>", press, add="+")
        canvas.bind_all("<B1-Motion>", motion, add="+")
        canvas.bind_all("<ButtonRelease-1>", release, add="+")
        canvas.bind_all("<MouseWheel>", lambda e: canvas.yview_scroll(int(-e.delta / 120), "units"), add="+")

    def _header(self):
        f = tk.Frame(self.frame, bg="#263238", padx=8, pady=6)
        f.pack(fill="x", padx=6, pady=(6, 5))
        bf = tk.Frame(f, bg="#263238")
        bf.pack(side="right", padx=(5, 0))
        tk.Button(bf, text="VERBETES TÉCNICOS", command=self.show_glossary,
                  font=("Arial", 8, "bold"), bg="#455a64", fg="white",
                  relief="flat", padx=5, pady=4).pack(fill="x")
        tk.Button(bf, text="VALIDAÇÃO 2019", command=self.show_validation,
                  font=("Arial", 8, "bold"), bg="#00695c", fg="white",
                  relief="flat", padx=5, pady=4).pack(fill="x", pady=(4, 0))

        title = tk.Frame(f, bg="#263238")
        title.pack(side="left", fill="x", expand=True)
        tk.Label(title, text="ZAB", fg="white", bg="#263238",
                 font=("Arial", 18, "bold"), anchor="w").pack(fill="x")
        tk.Label(title, text="ZnCl₂ / NH₄Cl / GDE  •  modelo v2 (quasipartículas, carga/descarga)",
                 fg="#e0e0e0", bg="#263238", anchor="w",
                 font=("Arial", 8)).pack(fill="x", pady=(1, 0))

    def _card(self, title, parent=None):
        box = tk.Frame(parent or self.frame, bg="#ffffff", bd=1, relief="solid")
        box.pack(fill="x", padx=6, pady=4)
        tk.Label(box, text=title, bg="#ffffff", anchor="w", justify="left", wraplength=330,
                 font=("Arial", 10, "bold")).pack(fill="x", padx=9, pady=(7, 4))
        return box

    def add_entry(self, parent, label, value):
        row = tk.Frame(parent, bg="#ffffff")
        row.pack(fill="x", padx=9, pady=2)
        tk.Label(row, text=label, width=20, anchor="w", bg="#ffffff",
                 font=("Arial", 9)).pack(side="left")
        v = tk.StringVar(value=str(value))
        self.vars[label] = v
        ent = tk.Entry(row, textvariable=v, width=12, justify="right",
                       font=("Arial", 9))
        ent.pack(side="right")

    def add_choice(self, parent, label, values, value):
        row = tk.Frame(parent, bg="#ffffff")
        row.pack(fill="x", padx=9, pady=2)
        tk.Label(row, text=label, width=20, anchor="w", bg="#ffffff",
                 font=("Arial", 9)).pack(side="left")
        v = tk.StringVar(value=str(value))
        self.vars[label] = v
        cb = ttk.Combobox(row, textvariable=v, values=list(values), width=22 if label == "Predefinição" else 10,
                          state="readonly", font=("Arial", 9))
        cb.pack(side="right")
        return cb

    def _input_card(self, card):
        box = self._card("CÉLULA E ELETRÓLITO")
        cb = self.add_choice(box, "Predefinição", list(PRESETS.keys()), "Eletrólito C — célula padrão")
        cb.bind("<<ComboboxSelected>>", self.apply_preset)
        self.add_entry(box, "ZnCl₂ (mol/L)", 0.5)
        self.add_entry(box, "NH₄Cl (mol/L)", 1.6)
        self.add_entry(box, "pH inicial", 8.0)
        self.add_choice(box, "Modelo de Zn", ["paste", "foil"], "paste")
        self.add_choice(box, "Sólidos", ["todos", "só ZnO"], "todos")
        self.add_entry(box, "Zn (mm)", 5.0)
        self.add_entry(box, "Separador (mm)", 0.5)
        self.add_entry(box, "GDE (mm)", 0.5)
        tk.Label(box, text="Zn (mm) só vale para 'paste'; 'foil' = folha, reação na superfície.",
                 bg="#ffffff", fg="#666666", font=("Arial", 8)).pack(anchor="w", padx=10, pady=(0, 7))

    def _protocol_card(self, card):
        box = self._card("PROTOCOLO")
        self.add_entry(box, "Ciclos", 1)
        self.add_entry(box, "I descarga (mA/cm²)", 1.0)
        self.add_entry(box, "t máx descarga (h)", 2000)
        self.add_entry(box, "V corte descarga (V)", 0.6)
        self.add_entry(box, "I carga (mA/cm²)", 0.0)
        self.add_entry(box, "t máx carga (h)", 8)
        self.add_entry(box, "V corte carga (V)", 2.6)
        self.add_entry(box, "Espessuras Zn (mm)", "0.5, 2, 5")
        tk.Label(box, text="I carga = 0: só descarga. Passo de tempo é adaptativo (sem 'passos × dt'). 'Espessuras Zn' só vale para a VARREDURA.",
                 bg="#ffffff", fg="#666666", font=("Arial", 8)).pack(anchor="w", padx=10, pady=(0, 7))

    def apply_preset(self, event=None):
        p = PRESETS.get(self.vars["Predefinição"].get())
        if not p:
            return
        for k, v in p.items():
            self.vars[k].set(str(v))

    def _run_card(self, card):
        box = self._card("EXECUÇÃO")
        # Barra própria em Canvas: mais visível e confiável no Tkinter/Android
        self.progress_canvas = tk.Canvas(
            box, height=34, bg="#e6e6e6", highlightthickness=1,
            highlightbackground="#bdbdbd", bd=0
        )
        self.progress_canvas.pack(fill="x", padx=10, pady=(3, 3))
        self.progress_fill = self.progress_canvas.create_rectangle(
            0, 0, 0, 34, fill="#1976d2", outline=""
        )
        self.progress_text = self.progress_canvas.create_text(
            0, 17, text="0.0%   |   segmento 0/0",
            fill="#222222", font=("Arial", 9, "bold"), anchor="center"
        )
        self.progress_canvas.bind("<Configure>", self._resize_progress)

        self.status = tk.StringVar(value="Pronto para iniciar")
        tk.Label(box, textvariable=self.status, bg="#ffffff", fg="#333333",
                 anchor="w", justify="left", font=("Arial", 8), wraplength=330).pack(
                     fill="x", padx=10, pady=(0, 4))

        buttons = tk.Frame(box, bg="#ffffff")
        buttons.pack(fill="x", padx=8, pady=(0, 8))
        self.start_btn = tk.Button(buttons, text="INICIAR", command=self.start,
                                   font=("Arial", 10, "bold"), height=1)
        self.start_btn.pack(side="left", fill="x", expand=True, padx=(2, 3))
        self.stop_btn = tk.Button(buttons, text="PARAR", command=self.stop,
                                  font=("Arial", 9), height=1, state="disabled")
        self.stop_btn.pack(side="left", fill="x", expand=True, padx=3)
        self.clear_btn = tk.Button(buttons, text="LIMPAR", command=self.clear,
                                   font=("Arial", 9), height=1)
        self.clear_btn.pack(side="left", fill="x", expand=True, padx=(3, 2))
        self._sweep_button(box)

    def _summary_card(self, card):
        box = self._card("RESULTADOS PRINCIPAIS")
        names = [
            ("tempo_sim", "Tempo de simulação"),
            ("tempo_cpu", "Tempo computacional"),
            ("voltage", "Tensão final"),
            ("capacity", "Capacidade (descargas)"),
            ("energy", "Energia (descargas)"),
            ("ph", "pH final (GDE)"),
            ("phclass", "Estabilidade de pH"),
            ("main", "Precipitado predominante"),
            ("prod", "Classificação do produto"),
        ]
        for key, label in names:
            row = tk.Frame(box, bg="#ffffff")
            row.pack(fill="x", padx=9, pady=1)
            tk.Label(row, text=label, bg="#ffffff", fg="#555555",
                     anchor="w", font=("Arial", 8)).pack(side="left", fill="x", expand=True)
            v = tk.StringVar(value="—")
            self.summary_vars[key] = v
            tk.Label(row, textvariable=v, bg="#ffffff", fg="#111111",
                     anchor="e", justify="right", font=("Arial", 8, "bold"),
                     wraplength=170).pack(side="right", padx=(4, 0))

        self.energy_note = tk.StringVar(value="Energia: —")
        tk.Label(box, textvariable=self.energy_note, bg="#ffffff", fg="#555555",
                 anchor="w", font=("Arial", 8)).pack(fill="x", padx=9, pady=(4, 7))

    def _result_card(self, card):
        box = self._card("RESULTADO FINAL — RESUMO")
        result_frame = tk.Frame(box, bg="#fafafa")
        result_frame.pack(fill="both", expand=True, padx=7, pady=(0, 4))
        self.out = tk.Text(result_frame, height=24, width=40, wrap="word",
                           font=("Courier", 8), bg="#fafafa", bd=0,
                           relief="flat", padx=6, pady=5,
                           yscrollcommand=lambda first, last: result_scroll.set(first, last))
        result_scroll = tk.Scrollbar(result_frame, orient="vertical", command=self.out.yview)
        self.out.pack(side="left", fill="both", expand=True)
        result_scroll.pack(side="right", fill="y")
        self._deterministic_report = ""
        tk.Label(box,
                 text="A tela mostra os indicadores essenciais. O TXT contém o relatório técnico completo, incluindo os perfis espaciais.",
                 bg="#f2f2f2", fg="#444444", justify="left", anchor="w",
                 wraplength=330, font=("Arial", 8)).pack(fill="x", padx=8, pady=(2, 5))
        tk.Button(box, text="GRAVAR RELATÓRIO TÉCNICO COMPLETO EM TXT",
                  command=self.save_deterministic_txt,
                  font=("Arial", 9, "bold"), height=2).pack(fill="x", padx=7, pady=(0, 7))

    # ---------- verbetes técnicos (portado do final5; conteúdo reescrito para o v2) ----------
    def show_glossary(self):
        # Página própria dentro da janela principal; não usa janela suspensa.
        if hasattr(self, "glossary_page") and self.glossary_page.winfo_exists():
            self.glossary_page.lift()
            return

        self.glossary_page = tk.Frame(self.root, bg="#f2f2f2")
        self.glossary_page.place(x=0, y=0, relwidth=1, relheight=1)

        outer = tk.Frame(self.glossary_page, bg="#f2f2f2")
        outer.pack(fill="both", expand=True)
        canvas = tk.Canvas(outer, bg="#f2f2f2", highlightthickness=0)
        canvas.pack(fill="both", expand=True)
        frame = tk.Frame(canvas, bg="#f2f2f2")
        wid = canvas.create_window((0, 0), window=frame, anchor="nw")
        frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>", lambda e: canvas.itemconfigure(wid, width=max(300, e.width)))

        top = tk.Frame(frame, bg="#263238")
        top.pack(fill="x", padx=6, pady=(6, 4))
        tk.Button(top, text="← VOLTAR", command=self.close_glossary,
                  font=("Arial", 9, "bold"), bg="#455a64", fg="white",
                  relief="flat", padx=8, pady=5).pack(side="left", padx=6, pady=6)
        tk.Label(top, text="VERBETES TÉCNICOS", bg="#263238", fg="white",
                 font=("Arial", 14, "bold"), anchor="w"
                 ).pack(side="left", padx=(4, 6), pady=6)
        tk.Label(frame, text="Significado das grandezas e regras efetivamente presentes no código v2.",
                 bg="#f2f2f2", fg="#444444", justify="left", anchor="w",
                 font=("Arial", 8), wraplength=330).pack(fill="x", padx=10, pady=(0, 5))

        terms = GLOSSARY_TERMS

        # Filtro por seção para navegação dos verbetes.
        control = tk.Frame(frame, bg="#ffffff", bd=1, relief="solid")
        control.pack(fill="x", padx=6, pady=(0, 5))
        tk.Label(control, text="SEÇÃO", bg="#ffffff", fg="#555555",
                 font=("Arial", 8, "bold")).pack(anchor="w", padx=8, pady=(6, 1))
        section_names = ["TODAS AS SEÇÕES"] + sorted(dict.fromkeys(x[0] for x in terms), key=lambda x: int(x.split(" — ", 1)[0]))
        section_var = tk.StringVar(value=section_names[0])
        section_menu = tk.OptionMenu(control, section_var, *section_names)
        section_menu.configure(font=("Arial", 8), anchor="w", width=38)
        section_menu.pack(fill="x", padx=8, pady=(0, 4))
        count_var = tk.StringVar(value=f"{len(terms)} verbetes")
        tk.Label(control, textvariable=count_var, bg="#ffffff", fg="#546e7a",
                 font=("Arial", 8), anchor="w").pack(fill="x", padx=8, pady=(0, 6))

        results_frame = tk.Frame(frame, bg="#f2f2f2")
        results_frame.pack(fill="x", expand=True)

        def render_terms(*_):
            for child in results_frame.winfo_children():
                child.destroy()
            selected = section_var.get()
            shown = 0
            for sec, name, meaning in terms:
                if selected != "TODAS AS SEÇÕES" and sec != selected:
                    continue
                card = tk.Frame(results_frame, bg="#ffffff", bd=1, relief="solid")
                card.pack(fill="x", padx=6, pady=2)
                tk.Label(card, text=name, bg="#ffffff", fg="#111111", anchor="w",
                         font=("Arial", 9, "bold"), wraplength=320).pack(fill="x", padx=8, pady=(5, 1))
                tk.Label(card, text=sec, bg="#ffffff", fg="#607d8b", anchor="w",
                         font=("Arial", 7, "bold"), wraplength=320).pack(fill="x", padx=8, pady=(0, 1))
                tk.Label(card, text=meaning, bg="#ffffff", fg="#444444", anchor="w",
                         justify="left", wraplength=320, font=("Arial", 8)).pack(fill="x", padx=8, pady=(0, 5))
                if sec == "12 — REFERÊNCIAS BIBLIOGRÁFICAS":
                    url = ("https://arxiv.org/abs/1902.04400" if name.startswith("REF-03") else
                           "https://chemistry-europe.onlinelibrary.wiley.com/doi/full/10.1002/cssc.201701468"
                           if name.startswith("REF-01") else
                           "https://chemistry-europe.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fcssc.201701468&file=cssc201701468-sup-0001-misc_information.pdf")
                    tk.Button(card, text="ABRIR NA INTERNET", command=lambda u=url: webbrowser.open(u),
                              font=("Arial", 8, "bold"), height=1).pack(fill="x", padx=8, pady=(0, 6))
                shown += 1
            count_var.set(f"{shown} de {len(terms)} verbetes")
            canvas.configure(scrollregion=canvas.bbox("all"))

        section_var.trace_add("write", render_terms)
        render_terms()

        # Texto responsivo: ocupa toda a largura útil em retrato e paisagem.
        def resize_glossary(event=None):
            w = max(260, canvas.winfo_width() - 32)
            for card in results_frame.winfo_children():
                card.configure(width=w)
                for sub in card.winfo_children():
                    sub.configure(wraplength=max(220, w - 24))
        canvas.bind("<Configure>", lambda e: (canvas.itemconfigure(wid, width=max(300, e.width)), resize_glossary(e)))
        self._bind_touch_scroll(canvas)
        canvas.after(50, resize_glossary)

    def close_glossary(self):
        if hasattr(self, "glossary_page") and self.glossary_page.winfo_exists():
            self.glossary_page.destroy()

    # ---------- progresso / gráfico ----------
    def _resize_progress(self, event=None):
        w = max(1, self.progress_canvas.winfo_width())
        h = max(1, self.progress_canvas.winfo_height())
        self.progress_canvas.coords(self.progress_fill, 0, 0, 0, h)
        self.progress_canvas.coords(self.progress_text, w / 2, h / 2)

    def _set_progress(self, p, seg, total):
        p = max(0.0, min(1.0, float(p)))
        w = max(1, self.progress_canvas.winfo_width())
        h = max(1, self.progress_canvas.winfo_height())
        self.progress_canvas.coords(self.progress_fill, 0, 0, w * p, h)
        self.progress_canvas.coords(self.progress_text, w / 2, h / 2)
        self.progress_canvas.itemconfigure(
            self.progress_text,
            text=f"{100*p:5.1f}%   |   segmento {seg}/{total}"
        )
        self.progress_canvas.update_idletasks()

    def update_progress(self, p, seg, total, V, elapsed, t_h, hist, snap=None):
        self._set_progress(p, seg, total)
        self._plot_hist = hist
        if snap is not None:
            self._snap = snap
        self._tick += 1
        if self._tick % 3 == 0:
            self._draw_plot()
        self.status.set(
            f"t: {t_h:.4g} h   •   V: {V:.4f} V   •   CPU: {elapsed:.1f} s"
        )

    # ---------- gráficos (Canvas puro; sem matplotlib) ----------
    CHART_DEFS = [
        ("vt",  "TENSÃO E pH DO GDE × TEMPO"),
        ("vq",  "TENSÃO × CAPACIDADE"),
        ("ph",  "pH × TEMPO: Zn, GDE, mínimo e máximo"),
        ("eta", "PERDAS DE TENSÃO: η Zn, η GDE e ôhmica"),
        ("phx", "PERFIL DE pH × POSIÇÃO (estado atual)"),
        ("cx",  "PERFIL DE Zn~, NH3~ E NH4~ × POSIÇÃO"),
        ("esp", "ESPÉCIES DE Zn: lado Zn × GDE"),
        ("sol", "SÓLIDOS × TEMPO"),
        ("zn",  "ZINCO: fração média e casca δ"),
        ("dt",  "PASSO DE TEMPO × TEMPO"),
        ("sw",  "VARREDURA: ENERGIA E CAPACIDADE × ESPESSURA DE Zn"),
    ]
    REGION_FILL = {"Zn": "#f5ede8", "sep": "#eef3f8", "GDE": "#fdf3e6"}

    def _charts_cards(self, card):
        self.charts = {}
        self._plot_hist = []
        self._snap = None
        self._sweep_rows = []
        self._tick = 0
        self._fonts = {"t": tkfont.Font(family="Arial", size=7),
                       "l": tkfont.Font(family="Arial", size=8),
                       "b": tkfont.Font(family="Arial", size=8, weight="bold")}
        # todas as folgas escalam com a altura de linha da fonte (não são pixels fixos)
        self._ls = max(14, self._fonts["l"].metrics("linespace"))
        for key, title in self.CHART_DEFS:
            box = self._card(title)
            box.pack_configure(pady=8)
            cv = tk.Canvas(box, height=260, bg="#fafafa", highlightthickness=1,
                           highlightbackground="#bdbdbd", bd=0)
            cv.pack(fill="x", padx=10, pady=(3, 12))
            cv.bind("<Configure>", lambda e: self._draw_charts())
            self.charts[key] = cv

    @staticmethod
    def _ticks(lo, hi, n=4):
        span = hi - lo
        if span <= 0:
            return [lo]
        step = 10 ** math.floor(math.log10(span / n))
        for m in (1, 2, 5, 10):
            if span / (step * m) <= n:
                step *= m
                break
        t0 = math.ceil(lo / step - 1e-9) * step
        out, v = [], t0
        while v <= hi + 1e-9 * span:
            out.append(0.0 if abs(v) < 1e-12 * max(1.0, abs(hi)) else v)
            v += step
        return out

    @staticmethod
    def _fmt_ticks(ticks):
        """Mesmo número de casas decimais em todas as marcas do eixo."""
        d = 0
        while d < 6 and any(abs(round(t, d) - t) > 1e-9 * max(1.0, abs(t)) for t in ticks):
            d += 1
        return [f"{t:.{d}f}" for t in ticks]

    @staticmethod
    def _ds(xs, ys, maxn=220):
        n = len(xs)
        if n <= maxn:
            return list(xs), list(ys)
        k = int(math.ceil(n / maxn))
        ix = list(range(0, n, k))
        if ix[-1] != n - 1:
            ix.append(n - 1)
        return [xs[i] for i in ix], [ys[i] for i in ix]

    def _tw(self, cv, text, font):
        """Largura REAL do texto, medida pelo próprio Canvas (bbox de um item de texto temporário)."""
        i = cv.create_text(-5000, -5000, text=text, anchor="nw", font=font)
        b = cv.bbox(i)
        cv.delete(i)
        return float(b[2] - b[0]) if b else float(font.measure(text))

    def _legend_pos(self, cv, items, w):
        """Posições dos itens da legenda; larguras medidas no Canvas, com folga de segurança."""
        ls = self._ls; f = self._fonts["l"]
        pad = max(10, 0.8 * ls)
        sq, tg = max(9, 0.55 * ls), max(8, 0.55 * ls)
        gap = max(30, 2.4 * ls)                       # espaço horizontal mínimo entre itens
        pitch = max(28, 2.1 * ls)                     # passo vertical entre linhas da legenda
        x, row, pos = pad, 0, []
        for col, lab in items:
            iw = sq + tg + self._tw(cv, lab, f) * 1.08 + 2
            if x > pad and x + iw > w - pad:
                x, row = pad, row + 1
            pos.append((x, row, iw))
            x += iw + gap
        return pos, (row + 1) * pitch + 0.7 * ls, (sq, tg, pitch)

    def _draw_legend(self, cv, items, pos, y0, dims):
        f = self._fonts["l"]
        sq, tg, pitch = dims
        for (col, lab), (x, row, iw) in zip(items, pos):
            yc = y0 + 0.35 * self._ls + row * pitch + pitch / 2
            cv.create_rectangle(x, yc - sq / 2, x + sq, yc + sq / 2, fill=col, outline=col)
            cv.create_text(x + sq + tg, yc, text=lab, anchor="w", fill="#333333", font=f)

    def _plot(self, cv, series, logy=False, regions=None, xrange=None, markers=False,
              xlab="", yl=None, yr=None, empty="sem dados"):
        """series: (xs, ys, cor, rótulo, eixo 'l'|'r'). yl/yr: (texto, cor) do rótulo de unidade."""
        ft, fl, fb = self._fonts["t"], self._fonts["l"], self._fonts["b"]
        ls = self._ls
        cv.delete("all")
        w = max(240, cv.winfo_width())
        plot_h = max(150, 8.5 * ls)
        series = [s for s in series if len(s[0]) >= 1]
        if not series:
            cv.configure(height=plot_h * 0.45 + 40)
            cv.create_text(w / 2, (plot_h * 0.45 + 40) / 2, text=empty, fill="#888888", font=fl)
            return
        has_r = any(s[4] == "r" for s in series)
        xs_all = [x for s in series for x in s[0]]
        x0, x1 = xrange if xrange else (min(xs_all), max(xs_all))
        if x1 - x0 < 1e-12:
            x1 = x0 + 1.0

        def rng(axis):
            vals = [y for s in series if s[4] == axis for y in s[1] if (y > 0 or not (logy and axis == "l"))]
            if not vals:
                return None
            if logy and axis == "l":
                vals = [math.log10(v) for v in vals]
            lo, hi = min(vals), max(vals)
            minspan = max(0.015 * max(abs(lo), abs(hi)), 1e-9)      # evita rótulos com casas decimais demais
            if hi - lo < minspan:
                mid = 0.5 * (hi + lo)
                lo, hi = mid - 0.5 * minspan, mid + 0.5 * minspan
            if hi - lo < 1e-9:
                hi, lo = hi + 0.5, lo - 0.5
            pad = 0.06 * (hi - lo)
            return lo - pad, hi + pad

        rl, rr = rng("l"), rng("r")
        tl = self._ticks(*rl, 5) if rl else []
        trr = self._ticks(*rr, 5) if rr else []
        lab_l = ([f"1e{t:.0f}" for t in tl] if logy else self._fmt_ticks(tl)) if tl else []
        lab_r = self._fmt_ticks(trr) if trr else []
        wl = max([self._tw(cv, s, ft) for s in lab_l] + [0])
        wr = max([self._tw(cv, s, ft) for s in lab_r] + [0])
        L = max(0.9 * ls + 22, wl + 0.8 * ls)
        R = (wr + 0.8 * ls) if has_r else 0.8 * ls
        items = [(s[2], s[3]) for s in series]
        lpos, lh, ldims = self._legend_pos(cv, items, w)
        unit_h = 1.7 * ls if (yl or yr) else 0
        reg_h = 1.5 * ls if regions else 0
        T_ = 0.4 * ls + lh + unit_h + reg_h + 0.5 * ls
        B = 0.5 * ls + 1.2 * ls + 0.5 * ls + (1.4 * ls if xlab else 0) + 0.4 * ls
        H = T_ + plot_h + B
        cv.configure(height=int(H))
        pw, ph_ = w - L - R, plot_h
        X = lambda v: L + (v - x0) / (x1 - x0) * pw

        def Y(v, r, lg):
            if lg:
                v = math.log10(max(v, 1e-30))
            return T_ + ph_ - (v - r[0]) / (r[1] - r[0]) * ph_

        # legenda e rótulos de unidade: faixas separadas, sem tocar nas curvas
        self._draw_legend(cv, items, lpos, 0.4 * ls, ldims)
        ytxt = 0.4 * ls + lh + 0.85 * ls
        if yl:
            cv.create_text(8, ytxt, text=yl[0], anchor="w", fill=yl[1], font=fb)
        if yr:
            cv.create_text(w - 8, ytxt, text=yr[0], anchor="e", fill=yr[1], font=fb)
        # regiões: sombreado + nome na faixa superior (omitido se não couber; ancorado na borda se estreita)
        if regions:
            drawn = []
            ry = T_ - 0.75 * ls
            for ri, (ra, rb, name) in enumerate(regions):
                xa, xb = X(max(ra, x0)), X(min(rb, x1))
                if xb - xa <= 0:
                    continue
                cv.create_rectangle(xa, T_, xb, T_ + ph_, fill=self.REGION_FILL.get(name, "#f5f5f5"), outline="")
                tw = self._tw(cv, name, ft) * 1.08
                if xb - xa >= tw + 10:
                    cv.create_text((xa + xb) / 2, ry, text=name, anchor="center", fill="#607d8b", font=ft)
                    drawn.append(((xa + xb) / 2 - tw / 2, (xa + xb) / 2 + tw / 2))
                elif ri == len(regions) - 1:
                    if all(xb - tw > d1 + 10 for (d0, d1) in drawn):
                        cv.create_text(xb, ry, text=name, anchor="e", fill="#607d8b", font=ft)
                        drawn.append((xb - tw, xb))
                elif ri == 0:
                    cv.create_text(xa, ry, text=name, anchor="w", fill="#607d8b", font=ft)
                    drawn.append((xa, xa + tw))
            for (ra, rb, name) in regions[:-1]:
                if x0 < rb < x1:
                    cv.create_line(X(rb), T_, X(rb), T_ + ph_, fill="#9e9e9e", dash=(3, 3))
        cv.create_rectangle(L, T_, w - R, T_ + ph_, outline="#bdbdbd")
        # eixo x
        nx = max(2, int(pw // (4.2 * ls)))
        tx_ticks = self._ticks(x0, x1, nx)
        for tx, lab in zip(tx_ticks, self._fmt_ticks(tx_ticks)):
            cv.create_line(X(tx), T_, X(tx), T_ + ph_, fill="#ececec")
            cv.create_line(X(tx), T_ + ph_, X(tx), T_ + ph_ + 4, fill="#757575")
            cv.create_text(X(tx), T_ + ph_ + 0.5 * ls, text=lab, anchor="n", fill="#444444", font=ft)
        if xlab:
            cv.create_text(L + pw / 2, H - 0.4 * ls, text=xlab, anchor="s", fill="#666666", font=ft)
        # eixos y
        lcols = [s[2] for s in series if s[4] == "l"]
        lcol = lcols[0] if len(lcols) == 1 else "#333333"
        for tv, lab in zip(tl, lab_l):
            yy = Y(10 ** tv if logy else tv, rl, logy)
            cv.create_line(L, yy, w - R, yy, fill="#ececec")
            cv.create_line(L - 4, yy, L, yy, fill="#757575")
            cv.create_text(L - 7, yy, text=lab, anchor="e", fill=lcol, font=ft)
        rcols = [s[2] for s in series if s[4] == "r"]
        rcol = rcols[0] if len(rcols) == 1 else "#333333"
        for tv, lab in zip(trr, lab_r):
            yy = Y(tv, rr, False)
            cv.create_line(w - R, yy, w - R + 4, yy, fill="#757575")
            cv.create_text(w - R + 7, yy, text=lab, anchor="w", fill=rcol, font=ft)
        # curvas
        for xs, ys, col, lab, ax in series:
            r = rl if ax == "l" else rr
            lg = logy and ax == "l"
            if r is None:
                continue
            xd, yd = self._ds(xs, ys)
            pts = [q for a, b in zip(xd, yd) if (b > 0 or not lg) for q in (X(a), Y(b, r, lg))]
            if len(pts) >= 4:
                cv.create_line(*pts, fill=col, width=2)
            if markers or len(pts) == 2:
                for i in range(0, len(pts), 2):
                    cv.create_oval(pts[i] - 3.5, pts[i + 1] - 3.5, pts[i] + 3.5, pts[i + 1] + 3.5, fill=col, outline=col)

    def _draw_bars(self, cv):
        cv.delete("all")
        snap = self._snap
        ls = self._ls
        w = max(240, cv.winfo_width())
        ft, fl = self._fonts["t"], self._fonts["l"]
        if not snap:
            cv.configure(height=int(max(150, 8.5 * ls) * 0.45 + 40))
            cv.create_text(w / 2, (max(150, 8.5 * ls) * 0.45 + 40) / 2, text="sem dados", fill="#888888", font=fl)
            return
        zz, zg = snap["z_zn"], snap["z_gde"]
        names = [k for k in zz if max(zz[k], zg[k]) >= 0.02]
        names.sort(key=lambda k: -max(zz[k], zg[k]))
        names = names[:7]
        items = [("#8d6e63", "lado Zn"), ("#1976d2", "GDE")]
        lpos, lh, ldims = self._legend_pos(cv, items, w)
        rowh = max(34, 2.4 * ls)
        bh = max(9, 0.6 * ls)
        H = 0.4 * ls + lh + max(1, len(names)) * rowh + 0.6 * ls
        cv.configure(height=int(H))
        self._draw_legend(cv, items, lpos, 0.4 * ls, ldims)
        y0 = 0.4 * ls + lh
        if not names:
            cv.create_text(w / 2, y0 + rowh / 2, text="nenhuma espécie ≥ 2% do Zn", fill="#888888", font=fl)
            return
        x0 = 12 + max(self._tw(cv, k, ft) for k in names) * 1.08 + 12
        vw = self._tw(cv, "0.00", ft) * 1.08
        bw = max(40, w - x0 - vw - 24)
        for i, k in enumerate(names):
            y = y0 + i * rowh
            cv.create_text(x0 - 10, y + rowh / 2, text=k, anchor="e", fill="#333333", font=ft)
            for j, (val, col) in enumerate(((zz[k], "#8d6e63"), (zg[k], "#1976d2"))):
                yy = y + (rowh - 2 * bh - 4) / 2 + j * (bh + 4)
                cv.create_rectangle(x0, yy, x0 + max(1, bw * val), yy + bh, fill=col, outline=col)
                if val >= 0.02:
                    cv.create_text(x0 + bw * val + 6, yy + bh / 2, text=f"{val:.2f}", anchor="w", fill="#555555", font=ft)

    def _draw_charts(self):
        if not hasattr(self, "charts"):
            return
        H = self._plot_hist
        snap = self._snap
        col = lambda k: [r.get(k, 0.0) for r in H]
        t = [r["t_h"] for r in H]
        sw = self._sweep_rows
        B_, O_, G_, R_, P_, K_, N_ = "#1976d2", "#ef6c00", "#2e7d32", "#c62828", "#6a1b9a", "#212121", "#8d6e63"
        sol = lambda i: [s[i] for s in col("solids_mol_m2")]
        jobs = {
            "vt": lambda c: self._plot(c, [(t, col("V"), B_, "V", "l"), (t, col("pH_gde"), O_, "pH GDE", "r")],
                                       xlab="t (h)", yl=("V", B_), yr=("pH", O_)),
            "vq": lambda c: self._plot(c, [(col("cap_mAh_cm2"), col("V"), B_, "V × Q", "l")],
                                       xlab="capacidade líquida (mAh/cm²)", yl=("V", B_)),
            "ph": lambda c: self._plot(c, [(t, col("pH_zn"), N_, "lado Zn", "l"), (t, col("pH_gde"), O_, "GDE", "l"),
                                           (t, col("pH_min"), "#90a4ae", "mínimo", "l"), (t, col("pH_max"), "#455a64", "máximo", "l")],
                                       xlab="t (h)", yl=("pH", "#333333")),
            "eta": lambda c: self._plot(c, [(t, col("eta_z"), N_, "η Zn", "l"), (t, col("ohm"), G_, "ôhmica", "l"),
                                            (t, col("eta_g"), O_, "η GDE", "r")],
                                        xlab="t (h)", yl=("V (η Zn, ôhmica)", "#333333"), yr=("V (η GDE)", O_)),
            "sol": lambda c: self._plot(c, [(t, sol(0), P_, "ZnO", "l"), (t, sol(1), R_, "simonkolleite", "l"),
                                            (t, sol(2), B_, "Zn(OH)2", "l"), (t, sol(3), K_, "ZnCl2·2NH3", "l")],
                                        xlab="t (h)", yl=("mol/m²", "#333333")),
            "zn": lambda c: self._plot(c, [(t, col("eps_zn_mean"), B_, "fração de Zn", "l"), (t, col("delta_um"), O_, "casca δ", "r")],
                                       xlab="t (h)", yl=("fração", B_), yr=("δ (µm)", O_)),
            "dt": lambda c: self._plot(c, [(t, col("dt"), "#455a64", "dt", "l")], logy=True,
                                       xlab="t (h)", yl=("dt (s)", "#455a64")),
            "sw": lambda c: self._plot(c, [([r["mm"] for r in sw], [r["Wh_L"] for r in sw], B_, "Wh/L", "l"),
                                           ([r["mm"] for r in sw], [r["cap"] for r in sw], O_, "mAh/cm²", "r")],
                                       markers=True, xlab="espessura de Zn (mm)", yl=("Wh/L", B_), yr=("mAh/cm²", O_)),
        }
        if snap:
            xs, ed = snap["x"], snap["edges"]
            nz, ns = snap["n_zn"], snap["n_sep"]
            e_zn = ed[nz - 1] if nz else 0.0
            e_sep = ed[nz + ns - 1]
            regs = ([("Zn", 0.0, e_zn)] if nz else []) + [("sep", e_zn, e_sep), ("GDE", e_sep, ed[-1])]
            regs = [(a, b, n) for (n, a, b) in regs]
            xr = (0.0, ed[-1])
            jobs["phx"] = lambda c: self._plot(c, [(xs, snap["pH"], O_, "pH", "l")], regions=regs, xrange=xr,
                                               xlab="posição x (mm)", yl=("pH", O_))
            jobs["cx"] = lambda c: self._plot(c, [(xs, snap["Zt"], B_, "Zn~", "l"), (xs, snap["A3"], G_, "NH3~", "l"),
                                                  (xs, snap["A4"], N_, "NH4~", "r")], regions=regs, xrange=xr,
                                              xlab="posição x (mm)", yl=("mol/L (Zn~, NH3~)", "#333333"), yr=("mol/L (NH4~)", N_))
        for key, cv in self.charts.items():
            try:
                if key == "esp":
                    self._draw_bars(cv)
                elif key == "zn" and getattr(self.model, "zn_model", "paste") == "foil":
                    self._plot(cv, [], empty="não se aplica (Zn em folha)")
                elif key == "sw" and not sw:
                    self._plot(cv, [], empty="use o botão VARREDURA")
                elif key in jobs:
                    jobs[key](cv)
                else:
                    self._plot(cv, [])
            except Exception as e:
                cv.delete("all")
                cv.create_text(6, 6, text=f"erro no gráfico: {e}", anchor="nw", fill="#c62828", font=self._fonts["l"])

    def _draw_plot(self):
        self._draw_charts()

    def btn_state(self, on):
        state = "normal" if on else "disabled"
        self.start_btn.configure(state=state)
        self.clear_btn.configure(state=state)
        self.sweep_btn.configure(state=state)
        self.stop_btn.configure(state="disabled" if on else "normal")

    def clear(self):
        if self.running:
            return
        self._set_progress(0.0, 0, 0)
        self.status.set("Pronto para iniciar")
        self.out.delete("1.0", "end")
        for v in self.summary_vars.values():
            v.set("—")
        self.energy_note.set("Energia: —")
        self._plot_hist = []
        self._snap = None
        self._sweep_rows = []
        self._draw_plot()
        self.canvas.yview_moveto(0)

    # ---------- execução ----------
    def _read_inputs(self):
        try:
            g = lambda k: float(self.vars[k].get())
            setup = {
                "zncl2_M": g("ZnCl₂ (mol/L)"), "nh4cl_M": g("NH₄Cl (mol/L)"), "initial_pH": g("pH inicial"),
                "zn_model": self.vars["Modelo de Zn"].get(), "zn_mm": g("Zn (mm)"),
                "separator_mm": g("Separador (mm)"), "gde_mm": g("GDE (mm)"),
            }
            if self.vars["Sólidos"].get() == "só ZnO":
                setup["solids"] = ["ZnO"]
            cycles = int(float(self.vars["Ciclos"].get()))
            idis, tdis, vdis = g("I descarga (mA/cm²)"), g("t máx descarga (h)"), g("V corte descarga (V)")
            ichg, tchg, vchg = g("I carga (mA/cm²)"), g("t máx carga (h)"), g("V corte carga (V)")
        except Exception:
            messagebox.showerror("Entrada", "Use valores numéricos válidos.")
            return None
        if cycles < 1 or idis <= 0 or tdis <= 0 or ichg < 0 or (ichg > 0 and tchg <= 0):
            messagebox.showerror("Entrada", "Ciclos ≥ 1, correntes ≥ 0 (descarga > 0) e tempos > 0.")
            return None
        if setup["zncl2_M"] <= 0 or setup["nh4cl_M"] <= 0 or setup["separator_mm"] <= 0 or setup["gde_mm"] <= 0:
            messagebox.showerror("Entrada", "Concentrações e espessuras devem ser positivas.")
            return None
        return setup, cycles, (idis, tdis, vdis), (ichg, tchg, vchg)

    def start(self):
        inp = self._read_inputs()
        if inp is None:
            return
        setup, cycles, (idis, tdis, vdis), (ichg, tchg, vchg) = inp
        segs = []
        for _ in range(cycles if ichg > 0 else 1):
            segs.append({"mode": "discharge", "I_mA_cm2": idis, "t_max_h": tdis, "V_cut": vdis})
            if ichg > 0:
                segs.append({"mode": "charge", "I_mA_cm2": ichg, "t_max_h": tchg, "V_cut": vchg})

        self._set_progress(0.0, 0, len(segs))
        self.out.delete("1.0", "end")
        self._plot_hist = []
        self._snap = None
        self._sweep_rows = []
        self._draw_plot()
        self.status.set("Iniciando simulação...")
        self.btn_state(False)
        self.running = True
        th = threading.Thread(target=self.worker, args=(setup, segs), daemon=True)
        th.start()

    # ---------- varredura da espessura do Zn (densidade de energia, Fig. 9 do artigo) ----------
    def _sweep_button(self, box):
        self.sweep_btn = tk.Button(box, text="VARREDURA: ESPESSURA DE Zn × ENERGIA", command=self.start_sweep,
                                   font=("Arial", 9), height=1)
        self.sweep_btn.pack(fill="x", padx=10, pady=(0, 8))

    def start_sweep(self):
        inp = self._read_inputs()
        if inp is None:
            return
        setup, cycles, (idis, tdis, vdis), _chg = inp
        try:
            ths = [float(x) for x in self.vars["Espessuras Zn (mm)"].get().replace(";", ",").split(",") if x.strip()]
            assert 1 <= len(ths) <= 8 and all(x > 0 for x in ths)
        except Exception:
            messagebox.showerror("Entrada", "Espessuras Zn: de 1 a 8 valores positivos separados por vírgula.")
            return
        self._set_progress(0.0, 0, len(ths))
        self.out.delete("1.0", "end")
        self._plot_hist = []
        self._snap = None
        self._sweep_rows = []
        self._draw_plot()
        self.status.set("Iniciando varredura...")
        self.btn_state(False)
        self.running = True
        th = threading.Thread(target=self._sweep_worker, args=(setup, ths, idis, tdis, vdis), daemon=True)
        th.start()

    def _sweep_worker(self, setup, ths, idis, tdis, vdis):
        rows = []
        try:
            for i, mm in enumerate(ths):
                st = dict(setup)
                st["zn_model"] = "paste"; st["zn_mm"] = mm; st["n_zn"] = int(max(16, 4*mm)); st["dt_max"] = 300.0
                model = Model(st)
                self.model = model

                def cb(m, k, n, frac, V, cpu, i=i):
                    self.root.after(0, self.update_progress, (i + frac)/len(ths), i + 1, len(ths), V, cpu, m.t/3600.0, [])

                r = model.run_protocol([{"mode": "discharge", "I_mA_cm2": idis, "t_max_h": tdis, "V_cut": vdis}],
                                       rec_every_s=1e9, callback=cb)[0]
                rows.append({"mm": mm, "Wh_L": model.energy_density_WhL(), "cap": r["cap_mAh_cm2"], "reason": r["reason"],
                             "pH": r["pH_gde"], "cpu": model.cpu_s, "V": model.V})
                self.root.after(0, self._sweep_partial, list(rows))
                if model.stop_flag:
                    break
            self.root.after(0, self._sweep_finish, rows, None)
        except Exception as e:
            self.root.after(0, self._sweep_finish, rows, e)

    def _sweep_partial(self, rows):
        self._sweep_rows = rows
        self._draw_plot()

    def _sweep_finish(self, rows, error):
        self.running = False
        t = ["VARREDURA — ESPESSURA DE Zn × DENSIDADE DE ENERGIA", "=" * 44,
             "Descarga até o corte; Wh por litro do volume total da célula.", "",
             "Zn(mm)  E(Wh/L)  Q(mAh/cm²)  pH GDE  CPU(s)  fim", "-" * 44]
        for r in rows:
            t.append(f"{r['mm']:6.3g} {r['Wh_L']:9.1f} {r['cap']:10.1f} {r['pH']:7.2f} {r['cpu']:7.1f}  {r['reason']}")
        if error is not None:
            t += ["", f"ERRO após {len(rows)} ponto(s): {error}"]
        self._deterministic_report = "\n".join(t)
        self.out.delete("1.0", "end")
        self.out.insert("1.0", self._deterministic_report)
        self._plot_hist = []
        self._snap = None
        self._sweep_rows = rows
        self._draw_plot()
        self.status.set(f"varredura concluída   •   {len(rows)} ponto(s)   •   (ver gráfico da varredura)")
        self._set_progress(1.0, len(rows), len(rows))
        self.canvas.yview_moveto(0)
        self.btn_state(True)

    def stop(self):
        if self.model is not None:
            self.model.stop_flag = True
            self.status.set("Interrompendo...")

    def worker(self, setup, segs):
        model = None
        try:
            model = Model(setup)
            self.model = model
            total_s = sum(sg["t_max_h"] for sg in segs)*3600.0

            def cb(m, k, n, frac, V, cpu):
                p = (k + frac)/n
                self.root.after(0, self.update_progress, p, min(k + 1, n), n, V, cpu, m.t/3600.0, list(m.history), m.snapshot())

            res = model.run_protocol(segs, rec_every_s=max(30.0, total_s/500.0), callback=cb)
            self.root.after(0, self.finish, model, res, None)
        except Exception as e:
            self.root.after(0, self.finish, model, None, e)

    # ---------- relatório ----------
    @staticmethod
    def _screen_summary(full_report):
        """Resumo legível em tela estreita; preserva o relatório integral para TXT."""
        cut = full_report.find("E. PERFIS ESPACIAIS")
        return full_report if cut < 0 else full_report[:cut].rstrip() + "\n\n(perfis espaciais: ver TXT)"

    def _build_report(self, model, res, elapsed):
        phclass, prodclass, eclass, main, ph, phspan, solid = classify(model)
        t = []
        t.append("RESULTADO FINAL — ZAB NEUTRA (v2)")
        t.append("=" * 44)
        t.append("REFERÊNCIA CIENTÍFICA: Clark / Latz / Horstmann — 2017")
        t.append("Implementação Python 1D simplificada (quasipartículas)")
        t.append("Não é reprodução integral do sistema DAE original em MATLAB")
        t.append("")
        t.append("A. IDENTIFICAÇÃO E EXECUÇÃO")
        t.append("-" * 44)
        t.append(f"Passos aceitos         : {model.n_steps}")
        t.append(f"Passos rejeitados      : {model.n_rej}")
        t.append(f"Tempo de simulação     : {model.t:.6g} s ({model.t/3600:.6g} h)")
        t.append(f"Tempo computacional    : {elapsed:.3f} s")
        t.append(f"Volumes espaciais      : {model.N}")
        t.append(f"Modelo de Zn           : {model.zn_model}")
        t.append(f"Separador / GDE        : {model.L_SEP*1000:.6g} / {model.L_GDE*1000:.6g} mm")
        t.append(f"pH inicial configurado : {model.initial_pH:.6f}")
        d = (model.invariants() - model.inv0)/model.inv0
        t.append(f"Deriva Zn/Cl/N         : {d[0]:.2e} / {d[1]:.2e} / {d[2]:.2e}")
        t.append("")
        t.append("B. DESEMPENHO ELETROQUÍMICO")
        t.append("-" * 44)
        t.append(f"Tensão final           : {model.V:.6f} V")
        t.append(f"pH final (GDE)         : {ph:.5f}   (mín/máx no ensaio: {model.pH_min_run:.3f}/{model.pH_max_run:.3f})")
        t.append(f"Estabilidade de pH     : {phclass}")
        t.append(f"Precipitado predominante: {main}")
        t.append(f"Classificação          : {prodclass}")
        t.append("")
        t.append("C. SEGMENTOS")
        t.append("-" * 44)
        for i, r in enumerate(res, 1):
            t.append(f"{i:>2} {r['mode']:<9} I={r['I_mA_cm2']:+.4g} mA/cm²  t={r['t_h']:.4g} h  Q={r['cap_mAh_cm2']:.4g} mAh/cm²")
            t.append(f"     E={r['E_Wh_m2']:.4g} Wh/m²  V_fim={r['V_end']:.4f} V  pH GDE={r['pH_gde']:.3f}  fim: {r['reason']}")
        t.append("")
        t.append("D. SÓLIDOS (mol/m² de célula)")
        t.append("-" * 44)
        for k in range(4):
            t.append(f"{SOLID_LABELS[k]:<22}: {solid[k]:.6g}")
        t.append("")
        t.append("E. PERFIS ESPACIAIS (estado final)")
        t.append("-" * 44)
        t.append("x(mm)   reg  pH     Zn~(M)    NH3~(M)   NH4~(M)   eps_Zn")
        ph_all = model.pH()
        for i in range(model.N):
            t.append(f"{model.x[i]*1000:7.3f} {model.region[i]:<4} {ph_all[i]:6.3f} {model.Q[i,0]/1000:9.4g} "
                     f"{model.Q[i,1]/1000:9.4g} {model.Q[i,2]/1000:9.4g} {model.eps_zn[i]:6.3f}")
        t.append("")
        t.append("F. PREMISSAS")
        t.append("-" * 44)
        t.append("Concentrações no lugar de atividades; log Ksp dos hidróxidos +%.1f (calibração);" % KSP_SHIFT_OH)
        t.append(f"i0_orr={model.i0_orr:.3g} A/m², i0_oer={model.i0_oer:.3g} A/m² (calibrados em curvas digitalizadas);")
        t.append("sem inundação do GDE, sem balanço de H2O, sem dendritos, pO2 constante.")
        return "\n".join(t), (phclass, prodclass, main, ph)

    def finish(self, model, res, error):
        self.running = False
        if error is not None:
            self.status.set("Erro na simulação")
            self.out.delete("1.0", "end")
            self.out.insert("1.0", f"ERRO:\n{error}\n")
            self.btn_state(True)
            return
        elapsed = getattr(model, "cpu_s", 0.0)
        report, (phclass, prodclass, main, ph) = self._build_report(model, res, elapsed)
        dis = [r for r in res if r["mode"] == "discharge"]
        cap = sum(r["cap_mAh_cm2"] for r in dis)*10.0            # mAh/cm² -> Ah/m²
        en = sum(r["E_Wh_m2"] for r in dis)
        self.summary_vars["tempo_sim"].set(f"{model.t:.6g} s")
        self.summary_vars["tempo_cpu"].set(f"{elapsed:.3f} s")
        self.summary_vars["voltage"].set(f"{model.V:.6f} V")
        self.summary_vars["capacity"].set(f"{cap:.6g} Ah/m²")
        self.summary_vars["energy"].set(f"{en:.6g} Wh/m²")
        self.summary_vars["ph"].set(f"{ph:.5f}")
        self.summary_vars["phclass"].set(phclass)
        self.summary_vars["main"].set(main)
        self.summary_vars["prod"].set(prodclass)
        if dis:
            last = dis[-1]
            self.energy_note.set(f"Última descarga: {last['E_Wh_m2']/model.dx.sum()/1000:.4g} Wh/L (volume total da célula). Fim: {res[-1]['reason']}.")
        self._deterministic_report = report
        self.out.delete("1.0", "end")
        self.out.insert("1.0", self._screen_summary(report))
        self.out.yview_moveto(0)
        self._plot_hist = list(model.history)
        self._snap = model.snapshot()
        self._sweep_rows = []
        self._draw_plot()
        self._set_progress(1.0, len(res), len(res))
        self.status.set(f"concluído   •   t={model.t/3600:.4g} h   •   V={model.V:.4f} V   •   fim: {res[-1]['reason']}")
        self.canvas.yview_moveto(0)
        self.btn_state(True)

    def save_deterministic_txt(self):
        report = getattr(self, "_deterministic_report", "")
        if not report:
            self.status.set("Nenhum resultado disponível para gravação.")
            return

        paths = [
            "/storage/emulated/0/Download/ZAB_resultado_v2.txt",
            "/sdcard/Download/ZAB_resultado_v2.txt",
        ]
        saved = None
        for path in paths:
            try:
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with open(path, "w", encoding="utf-8") as f:
                    f.write(report)
                saved = path
                break
            except Exception:
                pass

        if saved:
            self.status.set(f"TXT gravado: {saved}")
        else:
            self.status.set("TXT não pôde ser gravado no armazenamento Android.")


    # ============================================================
    # ETAPA 2 — VALIDAÇÃO EXPERIMENTAL (2019): página própria
    # ============================================================
    VAL_CHART_DEFS = [
        ("k",   "CONDUTIVIDADE IÔNICA: modelo × artigo (mS/cm)", "execute o teste 1"),
        ("orr", "TENSÃO DE DESCARGA (ORR): modelo × artigo (V)", "execute o teste 2"),
        ("oer", "TENSÃO DE CARGA (OER): modelo × artigo (V)", "execute o teste 2"),
        ("vt",  "TENSÃO × TEMPO no ciclo 1 (1 mA/cm², 2 h + 2 h)", "execute o teste 2"),
        ("par", "pH DO AR × TEMPO (operando, 2 mA/cm²)", "execute o teste 3"),
        ("pzn", "pH DO LADO Zn × TEMPO (operando, 2 mA/cm²)", "execute o teste 3"),
        ("fas", "FASES SÓLIDAS APÓS 40% DE DoD: modelo × artigo", "execute o teste 4"),
        ("znb", "BALANÇO DE Zn: dissolvido e precipitado (E6)", "execute o teste 5"),
    ]
    VAL_NOTA_XRD = {"E4": "simonkolleíta", "E6": "simonkolleíta", "E7": "ZnCl₂·2NH₃", "E8": "camada intercalada (não modelada)"}

    def _stages_card(self):
        box = self._card("ETAPAS DO TRABALHO")
        rows = (("1", "SIMULADOR — modelo v2", "Reconstrução do modelo de 2017 (Clark, Latz, Horstmann): esta página."),
                ("2", "VALIDAÇÃO EXPERIMENTAL — 2019", "O mesmo modelo, sem recalibração, contra os experimentos de Clark, Mainar et al."))
        for n, tit, desc in rows:
            r = tk.Frame(box, bg="#ffffff")
            r.pack(fill="x", padx=9, pady=2)
            tk.Label(r, text=n, bg="#263238", fg="white", width=2, font=("Arial", 10, "bold")).pack(side="left", anchor="n", padx=(0, 7))
            col = tk.Frame(r, bg="#ffffff")
            col.pack(side="left", fill="x", expand=True)
            tk.Label(col, text=tit, bg="#ffffff", anchor="w", font=("Arial", 8, "bold")).pack(fill="x")
            tk.Label(col, text=desc, bg="#ffffff", fg="#555555", anchor="w", justify="left", wraplength=270,
                     font=("Arial", 8)).pack(fill="x")
        tk.Button(box, text="ABRIR ETAPA 2 — VALIDAÇÃO 2019", command=self.show_validation,
                  font=("Arial", 9, "bold"), height=1).pack(fill="x", padx=9, pady=(4, 8))

    def _val_alive(self):
        try:
            return hasattr(self, "val_page") and bool(self.val_page.winfo_exists())
        except Exception:
            return False

    def show_validation(self):
        if self._val_alive():
            self.val_page.lift()
            return
        self.val_page = tk.Frame(self.root, bg="#f2f2f2")
        self.val_page.place(x=0, y=0, relwidth=1, relheight=1)
        outer = tk.Frame(self.val_page, bg="#f2f2f2")
        outer.pack(fill="both", expand=True)
        canvas = tk.Canvas(outer, bg="#f2f2f2", highlightthickness=0)
        canvas.pack(fill="both", expand=True)
        frame = tk.Frame(canvas, bg="#f2f2f2")
        wid = canvas.create_window((0, 0), window=frame, anchor="nw")
        frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>", lambda e: canvas.itemconfigure(wid, width=max(300, e.width)))
        self._vframe = frame
        top = tk.Frame(frame, bg="#263238")
        top.pack(fill="x", padx=6, pady=(6, 4))
        tk.Button(top, text="← VOLTAR", command=self.close_validation, font=("Arial", 9, "bold"),
                  bg="#455a64", fg="white", relief="flat", padx=8, pady=5).pack(side="left", padx=6, pady=6)
        tk.Label(top, text="ETAPA 2 — VALIDAÇÃO 2019", bg="#263238", fg="white", font=("Arial", 13, "bold"),
                 anchor="w").pack(side="left", padx=(4, 6), pady=6)
        self._val_context_card()
        self._val_run_card()
        self._val_board_card()
        self._val_charts_cards()
        self._val_assume_card()
        self._val_report_card()
        self._bind_touch_scroll(canvas)
        self._val_render()

    def close_validation(self):
        if self._val_alive():
            self.val_page.destroy()

    def _vcard(self, title):
        return self._card(title, parent=self._vframe)

    def _val_context_card(self):
        box = self._vcard("CONTEXTO")
        txt = ("O modelo de 2017 (Clark, Latz, Horstmann) previu pH ácido durante a carga e produtos de descarga diferentes de ZnO. "
               "Em 2019, Clark, Mainar et al. mediram isso em células reais (arXiv:1902.04400). Esta etapa roda no modelo v2 as mesmas "
               "condições e compara com o artigo, sem recalibrar nada: os desvios são o resultado.\n\n"
               "Eletrólitos do artigo: E6 = A, E7 = B, E8 = C; o E4 (pH 4) é novo.\n"
               "Células: C-EI (eletrodos a 0,9 cm) e C-OpH (2,8 cm); folha de Zn de 250 µm, área de 1,327 cm².")
        tk.Label(box, text=txt, bg="#ffffff", fg="#333333", anchor="w", justify="left", wraplength=320,
                 font=("Arial", 8)).pack(fill="x", padx=9, pady=(0, 8))

    def _val_run_card(self):
        box = self._vcard("EXECUÇÃO")
        st = self.val_state
        self.val_sel = {}
        for tid, title in VAL_TESTS:
            v = tk.BooleanVar(value=True)
            self.val_sel[tid] = v
            tk.Checkbutton(box, text=title, variable=v, bg="#ffffff", anchor="w", justify="left",
                           font=("Arial", 8), wraplength=310).pack(fill="x", padx=9)
        row = tk.Frame(box, bg="#ffffff")
        row.pack(fill="x", padx=9, pady=(5, 2))
        tk.Label(row, text="Modo", width=20, anchor="w", bg="#ffffff", font=("Arial", 9)).pack(side="left")
        self.val_mode = tk.StringVar(value=st["modos"].get("t3", "reduzido"))
        ttk.Combobox(row, textvariable=self.val_mode, values=["completo", "reduzido"], width=10,
                     state="readonly", font=("Arial", 9)).pack(side="right")
        tk.Label(box, text=("Completo: condições do artigo (29 h + 29 h, 40% de DoD, 10 h + 10 h); no PC leva cerca de 10 min e no celular bem mais. "
                            "Reduzido: ciclos mais curtos, marcados INDICATIVO e sem valor de comparação. "
                            "Os testes 1 e 2 são sempre completos."),
                 bg="#ffffff", fg="#666666", anchor="w", justify="left", wraplength=320,
                 font=("Arial", 8)).pack(fill="x", padx=9, pady=(0, 4))
        self.val_pcv = tk.Canvas(box, height=34, bg="#e6e6e6", highlightthickness=1, highlightbackground="#bdbdbd", bd=0)
        self.val_pcv.pack(fill="x", padx=10, pady=(3, 3))
        self.val_pfill = self.val_pcv.create_rectangle(0, 0, 0, 34, fill="#1976d2", outline="")
        self.val_ptext = self.val_pcv.create_text(0, 17, text="0.0%", fill="#222222", font=("Arial", 9, "bold"), anchor="center")
        self.val_pcv.bind("<Configure>", lambda e: self._val_set_progress(self._val_p))
        self._val_p = 0.0
        self.val_status = tk.StringVar(value="Pronto para iniciar")
        tk.Label(box, textvariable=self.val_status, bg="#ffffff", fg="#333333", anchor="w", justify="left",
                 font=("Arial", 8), wraplength=330).pack(fill="x", padx=10, pady=(0, 4))
        bt = tk.Frame(box, bg="#ffffff")
        bt.pack(fill="x", padx=8, pady=(0, 8))
        self.val_start_btn = tk.Button(bt, text="INICIAR", command=self.val_start, font=("Arial", 10, "bold"), height=1)
        self.val_start_btn.pack(side="left", fill="x", expand=True, padx=(2, 3))
        self.val_stop_btn = tk.Button(bt, text="PARAR", command=self.val_stop, font=("Arial", 9), height=1, state="disabled")
        self.val_stop_btn.pack(side="left", fill="x", expand=True, padx=3)
        self.val_clear_btn = tk.Button(bt, text="LIMPAR", command=self.val_clear, font=("Arial", 9), height=1)
        self.val_clear_btn.pack(side="left", fill="x", expand=True, padx=(3, 2))

    def _val_board_card(self):
        box = self._vcard("PLACAR DE VALIDAÇÃO")
        self.val_sum = tk.StringVar(value="Nenhum teste executado.")
        tk.Label(box, textvariable=self.val_sum, bg="#ffffff", fg="#111111", anchor="w", justify="left",
                 wraplength=320, font=("Arial", 9, "bold")).pack(fill="x", padx=9, pady=(0, 4))
        self.val_board = tk.Frame(box, bg="#ffffff")
        self.val_board.pack(fill="x", padx=6, pady=(0, 8))

    def _val_charts_cards(self):
        self.val_charts = {}
        for key, title, _msg in self.VAL_CHART_DEFS:
            box = self._vcard(title)
            box.pack_configure(pady=8)
            cv = tk.Canvas(box, height=260, bg="#fafafa", highlightthickness=1, highlightbackground="#bdbdbd", bd=0)
            cv.pack(fill="x", padx=10, pady=(3, 12))
            cv.bind("<Configure>", lambda e: self._val_draw_charts())
            self.val_charts[key] = cv

    def _val_assume_card(self):
        box = self._vcard("PREMISSAS E CRITÉRIOS")
        a = VAL_ASSUME
        txt = (f"• Sem recalibração: ajustar nos mesmos dados e depois validar seria circular.\n"
               f"• Tensão: o valor do artigo deve cair na faixa do modelo ± {a['tol_V']:.2f} V. A Tabela 3 não diz se é média, início ou fim do semiciclo; "
               f"por isso o placar também mostra o sobrepotencial médio.\n"
               f"• Condutividade: razão modelo/artigo dentro de ±{100*a['tol_kappa']:.0f}%.\n"
               f"• pH do E4 ao fim da descarga entre {a['plato_pH'][0]} e {a['plato_pH'][1]} (artigo: perto de pH 6); queda ácida e salto: limiar {a['limiar_pH']:.1f}, o mesmo de INSTÁVEL.\n"
               f"• Fases: sem XRD quantitativo; usa-se a fração molar de Zn em cada sólido. 40% de DoD = {a['dod40_mAh_cm2']:.0f} mAh/cm²; célula C-EI assumida.\n"
               f"• Fig. 12: condições não informadas no texto; usa-se E6, C-EI, 1 mA/cm², 10 h + 10 h.\n"
               f"• Balanço de Zn: erro até {100*a['faraday_tol']:.0f}% do Zn ciclado.\n"
               f"• Fora do alcance: evaporação, inundação do GDE, dissolução do MnO₂ em meio ácido, HER e a degradação em 150 ciclos.")
        tk.Label(box, text=txt, bg="#ffffff", fg="#333333", anchor="w", justify="left", wraplength=320,
                 font=("Arial", 8)).pack(fill="x", padx=9, pady=(0, 8))

    def _val_report_card(self):
        box = self._vcard("RELATÓRIO DA VALIDAÇÃO")
        fr = tk.Frame(box, bg="#fafafa")
        fr.pack(fill="both", expand=True, padx=7, pady=(0, 4))
        self.val_out = tk.Text(fr, height=22, width=40, wrap="word", font=("Courier", 8), bg="#fafafa", bd=0,
                               relief="flat", padx=6, pady=5, yscrollcommand=lambda a, b: sc.set(a, b))
        sc = tk.Scrollbar(fr, orient="vertical", command=self.val_out.yview)
        self.val_out.pack(side="left", fill="both", expand=True)
        sc.pack(side="right", fill="y")
        tk.Button(box, text="GRAVAR RELATÓRIO DA VALIDAÇÃO EM TXT", command=self.val_save_txt,
                  font=("Arial", 9, "bold"), height=2).pack(fill="x", padx=7, pady=(4, 7))

    # ---------- execução ----------
    def _val_buttons(self, idle):
        if not self._val_alive():
            return
        self.val_start_btn.configure(state="normal" if idle else "disabled")
        self.val_clear_btn.configure(state="normal" if idle else "disabled")
        self.val_stop_btn.configure(state="disabled" if idle else "normal")

    def _val_set_progress(self, p, text=None):
        if not self._val_alive():
            return
        self._val_p = max(0.0, min(1.0, float(p)))
        w = max(1, self.val_pcv.winfo_width()); h = max(1, self.val_pcv.winfo_height())
        self.val_pcv.coords(self.val_pfill, 0, 0, w*self._val_p, h)
        self.val_pcv.coords(self.val_ptext, w/2, h/2)
        self.val_pcv.itemconfigure(self.val_ptext, text=text if text is not None else f"{100*self._val_p:5.1f}%")

    def val_start(self):
        st = self.val_state
        if st["running"]:
            return
        tests = [tid for tid, _t in VAL_TESTS if self.val_sel[tid].get()]
        if not tests:
            self.val_status.set("Selecione ao menos um teste.")
            return
        modo = self.val_mode.get()
        for tid in tests:
            st["res"].pop(tid, None)
            st["checks"] = [c for c in st["checks"] if c["t"] != tid]
            st["modos"][tid] = "completo" if tid in ("t1", "t2") else modo
        st["running"] = True
        st["tests"] = tests
        st["report"] = ""
        self.val_ctl = ValCtl()
        self._val_set_progress(0.0)
        self.val_status.set("Iniciando validação...")
        self._val_buttons(False)
        self._val_render()
        threading.Thread(target=self._val_worker, args=(tests, modo, self.val_ctl), daemon=True).start()

    def val_stop(self):
        if self.val_state["running"] and getattr(self, "val_ctl", None) is not None:
            self.val_ctl.stop = True
            self.val_status.set("Interrompendo...")

    def val_clear(self):
        st = self.val_state
        if st["running"]:
            return
        st.update({"res": {}, "checks": [], "modos": {}, "report": "", "tests": []})
        self._val_set_progress(0.0)
        if self._val_alive():
            self.val_status.set("Pronto para iniciar")
        self._val_render()

    def _val_worker(self, tests, modo, ctl):
        t0 = time.time(); err = None
        after = self.root.after
        ctl.on_progress = lambda p, lab, V: after(0, self._val_progress, p, lab, V)
        ctl.on_partial = lambda tid, r: after(0, self._val_partial, tid, dict(r))
        ctl.total_h = sum(val_hours(t, modo) for t in tests)
        ctl.done_h = 0.0
        try:
            for tid in tests:
                if ctl.stop:
                    break
                out, res = [], {}
                VAL_FUNCS[tid](ctl, modo if tid not in ("t1", "t2") else "completo", out, res)
                after(0, self._val_test_done, tid, res.get(tid), out)
        except Exception as e:
            err = e
        after(0, self._val_finish, err, time.time() - t0)

    def _val_progress(self, p, label, V):
        if not self._val_alive():
            return
        self._val_set_progress(p)
        tail = f"  •  V {V:.3f} V" if (V is not None and V == V) else ""
        self.val_status.set(f"{label}{tail}")

    def _val_partial(self, tid, r):
        self.val_state["res"][tid] = r
        if self._val_alive():
            self._val_draw_charts()

    def _val_test_done(self, tid, r, checks):
        st = self.val_state
        if r is not None:
            st["res"][tid] = r
        st["checks"] = [c for c in st["checks"] if c["t"] != tid] + checks
        self._val_render()

    def _val_finish(self, err, elapsed):
        st = self.val_state
        st["running"] = False
        st["elapsed"] = elapsed
        for tid in st["tests"]:
            if not any(c["t"] == tid for c in st["checks"]):
                st["checks"].append({"t": tid, "txt": "teste não executado", "st": "PENDENTE", "det": "validação interrompida"})
        if err is not None:
            st["checks"].append({"t": st["tests"][0] if st["tests"] else "t1", "txt": "falha numérica", "st": "PENDENTE",
                                 "det": f"{type(err).__name__}: {err}"})
        st["report"] = val_report(st)
        if self._val_alive():
            self._val_buttons(True)
            self._val_set_progress(1.0 if err is None and not (getattr(self, "val_ctl", None) and self.val_ctl.stop) else self._val_p)
            self.val_status.set(f"Concluído em {elapsed:.0f} s." if err is None else f"Erro: {err}")
        self._val_render()

    # ---------- exibição ----------
    def _val_render(self):
        if not self._val_alive():
            return
        st = self.val_state
        for ch in self.val_board.winfo_children():
            ch.destroy()
        cc = val_counts(st["checks"])
        self.val_sum.set("  ".join(f"{k} {v}" for k, v in cc.items() if v) or "Nenhum teste executado.")
        for tid, title in VAL_TESTS:
            xs = [x for x in st["checks"] if x["t"] == tid]
            if not xs:
                continue
            tk.Label(self.val_board, text=title, bg="#eceff1", fg="#263238", anchor="w", justify="left", wraplength=320,
                     font=("Arial", 8, "bold")).pack(fill="x", padx=2, pady=(6, 2))
            for x in xs:
                row = tk.Frame(self.val_board, bg="#ffffff")
                row.pack(fill="x", padx=2, pady=2)
                tk.Label(row, text=x["st"], bg=VAL_STATUS_COR.get(x["st"], "#9e9e9e"), fg="#ffffff", width=10,
                         font=("Arial", 7, "bold")).pack(side="left", anchor="n", padx=(0, 6))
                col = tk.Frame(row, bg="#ffffff")
                col.pack(side="left", fill="x", expand=True)
                tk.Label(col, text=x["txt"], bg="#ffffff", fg="#111111", anchor="w", justify="left", wraplength=235,
                         font=("Arial", 8)).pack(fill="x")
                if x["det"]:
                    tk.Label(col, text=x["det"], bg="#ffffff", fg="#666666", anchor="w", justify="left", wraplength=235,
                             font=("Arial", 7)).pack(fill="x")
        self.val_out.delete("1.0", "end")
        self.val_out.insert("1.0", st["report"] or "O relatório aparece ao fim da execução.")
        self._val_draw_charts()

    def _val_draw_charts(self):
        if not self._val_alive() or not hasattr(self, "val_charts"):
            return
        res = self.val_state["res"]
        jobs = {}
        r1, r2, r3, r4, r5 = (res.get(k) for k in ("t1", "t2", "t3", "t4", "t5"))
        if r1:
            els = list(r1)
            jobs["k"] = lambda c: self._val_pairs(c, [VAL_LABEL[e] for e in els], [r1[e][0] for e in els],
                                                  [r1[e][1] for e in els], "modelo", "artigo", "%.0f")
        if r2:
            els = list(r2)
            cats = [VAL_LABEL[e] for e in els]
            jobs["orr"] = lambda c: self._val_pairs(c, cats, [r2[e]["D"]["med"] for e in els], [VAL_TAB3[e][0] for e in els],
                                                    "modelo (média)", "artigo", "%.3f",
                                                    whisk={VAL_LABEL[e]: (r2[e]["D"]["vmin"], r2[e]["D"]["vmax"]) for e in els})
            jobs["oer"] = lambda c: self._val_pairs(c, cats, [r2[e]["C"]["med"] for e in els], [VAL_TAB3[e][1] for e in els],
                                                    "modelo (média)", "artigo", "%.3f",
                                                    whisk={VAL_LABEL[e]: (r2[e]["C"]["vmin"], r2[e]["C"]["vmax"]) for e in els})
            jobs["vt"] = lambda c: self._plot(c, [([h[0] for h in r2[e]["hist"]], [h[1] for h in r2[e]["hist"]],
                                                   VAL_COR[e], VAL_LABEL[e], "l") for e in els],
                                              xlab="t (h)", yl=("V", "#333333"))
        if r3:
            els = list(r3)
            jobs["par"] = lambda c: self._plot(c, [([h[0] for h in r3[e]["hist"]], [h[1] for h in r3[e]["hist"]],
                                                    VAL_COR[e], VAL_LABEL[e], "l") for e in els],
                                               xlab="t (h)", yl=("pH", "#333333"))
            jobs["pzn"] = lambda c: self._plot(c, [([h[0] for h in r3[e]["hist"]], [h[2] for h in r3[e]["hist"]],
                                                    VAL_COR[e], VAL_LABEL[e], "l") for e in els],
                                               xlab="t (h)", yl=("pH", "#333333"))
        if r4:
            els = list(r4)
            jobs["fas"] = lambda c: self._val_stack(c, [VAL_LABEL[e] for e in els], [r4[e]["f"] for e in els],
                                                    ["artigo: " + self.VAL_NOTA_XRD[e] for e in els])
        if r5 and r5.get("samples"):
            s5 = r5["samples"]
            jobs["znb"] = lambda c: self._plot(c, [([s[0] for s in s5], [s[1] for s in s5], "#1976d2", "Zn dissolvido", "l"),
                                                   ([s[0] for s in s5], [s[2] for s in s5], "#c62828", "Zn em sólidos", "r")],
                                               xlab="t (h)", yl=("mol/m²", "#1976d2"), yr=("mol/m²", "#c62828"))
        for key, title, msg in self.VAL_CHART_DEFS:
            cv = self.val_charts[key]
            try:
                if key in jobs:
                    jobs[key](cv)
                else:
                    self._plot(cv, [], empty=msg)
            except Exception as e:
                cv.delete("all")
                cv.create_text(6, 6, text=f"erro no gráfico: {e}", anchor="nw", fill="#c62828", font=self._fonts["l"])

    def _val_pairs(self, cv, cats, A, B, la, lb, fmt, whisk=None):
        """Barras horizontais pareadas (modelo × artigo); 'whisk' marca a faixa mín–máx do modelo."""
        cv.delete("all")
        ls = self._ls; w = max(240, cv.winfo_width()); ft = self._fonts["t"]
        items = [("#1976d2", la), ("#8d6e63", lb)]
        lpos, lh, ldims = self._legend_pos(cv, items, w)
        rowh = max(38, 2.6*ls); bh = max(9, 0.6*ls)
        H = 0.4*ls + lh + len(cats)*rowh + 0.6*ls
        cv.configure(height=int(H))
        self._draw_legend(cv, items, lpos, 0.4*ls, ldims)
        y0 = 0.4*ls + lh
        x0 = 12 + max(self._tw(cv, k, ft) for k in cats)*1.08 + 12
        vw = self._tw(cv, fmt % 8888.888, ft)*1.08
        tops = [v for v in list(A) + list(B) if v is not None] + [b for k in (whisk or {}) for b in whisk[k]]
        sc = max(40.0, w - x0 - vw - 24)/max(max(tops) if tops else 1.0, 1e-9)
        for i, k in enumerate(cats):
            y = y0 + i*rowh
            cv.create_text(x0 - 10, y + rowh/2, text=k, anchor="e", fill="#333333", font=ft)
            for j, (val, col) in enumerate(((A[i], "#1976d2"), (B[i], "#8d6e63"))):
                if val is None:
                    continue
                yy = y + (rowh - 2*bh - 4)/2 + j*(bh + 4)
                end = x0 + val*sc
                cv.create_rectangle(x0, yy, max(x0 + 1, end), yy + bh, fill=col, outline=col)
                if j == 0 and whisk and k in whisk:
                    lo_, hi_ = whisk[k]
                    yc = yy + bh/2
                    cv.create_line(x0 + lo_*sc, yc, x0 + hi_*sc, yc, fill="#111111", width=2)
                    cv.create_line(x0 + lo_*sc, yc - bh/2, x0 + lo_*sc, yc + bh/2, fill="#111111", width=2)
                    cv.create_line(x0 + hi_*sc, yc - bh/2, x0 + hi_*sc, yc + bh/2, fill="#111111", width=2)
                    end = max(end, x0 + hi_*sc)
                cv.create_text(end + 6, yy + bh/2, text=fmt % val, anchor="w", fill="#555555", font=ft)

    def _val_stack(self, cv, cats, fr, notes):
        """Barras empilhadas: fração molar de Zn por sólido; abaixo de cada barra, a fase observada no artigo."""
        cv.delete("all")
        ls = self._ls; w = max(240, cv.winfo_width()); ft = self._fonts["t"]
        items = [(VAL_FASE_COR[i], VAL_FASE_NOME[SOLIDS[i]]) for i in range(4)]
        lpos, lh, ldims = self._legend_pos(cv, items, w)
        bh = max(12, 0.9*ls); rowh = 2.6*ls + bh/2 + 18
        H = 0.4*ls + lh + len(cats)*rowh + 0.6*ls
        cv.configure(height=int(H))
        self._draw_legend(cv, items, lpos, 0.4*ls, ldims)
        y0 = 0.4*ls + lh
        x0 = 12 + max(self._tw(cv, k, ft) for k in cats)*1.08 + 12
        bw = max(40.0, w - x0 - 16)
        for i, k in enumerate(cats):
            y = y0 + i*rowh + 6
            cv.create_text(x0 - 10, y + bh/2, text=k, anchor="e", fill="#333333", font=ft)
            x = x0
            for j, v in enumerate(fr[i]):
                ww = bw*v
                if ww > 0.5:
                    cv.create_rectangle(x, y, x + ww, y + bh, fill=VAL_FASE_COR[j], outline="#ffffff")
                x += ww
            cv.create_text(x0, y + bh + 0.5*ls, text=notes[i], anchor="nw", fill="#555555", font=ft)

    def val_save_txt(self):
        rep = self.val_state.get("report", "")
        if not rep:
            if self._val_alive():
                self.val_status.set("Nenhum resultado disponível para gravação.")
            return
        paths = ["/storage/emulated/0/Download/ZAB_validacao_2019.txt", "/sdcard/Download/ZAB_validacao_2019.txt",
                 "ZAB_validacao_2019.txt"]
        saved = None
        for path in paths:
            try:
                d = os.path.dirname(path)
                if d:
                    os.makedirs(d, exist_ok=True)
                with open(path, "w", encoding="utf-8") as f:
                    f.write(rep)
                saved = path
                break
            except Exception:
                pass
        if self._val_alive():
            self.val_status.set(f"TXT gravado: {saved}" if saved else "TXT não pôde ser gravado.")


if __name__ == "__main__":
    root = tk.Tk()
    App(root)
    root.mainloop()
