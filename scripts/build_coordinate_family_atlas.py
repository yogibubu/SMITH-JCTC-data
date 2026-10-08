"""Build a clean vector atlas of the coordinate families used by the manuscript."""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle, FancyArrowPatch, FancyBboxPatch, Polygon
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch

OUT = Path(__file__).resolve().parents[1] / "figures"
NAVY = "#17324d"
TEXT = "#506274"
LINE = "#8395a2"
TEAL = "#15958f"
AMBER = "#d58a24"
VIOLET = "#7657a5"
BLUE = "#377bb4"
RED = "#c75a4a"
PALE = "#f7fafb"


def atom(ax, xy, color="#3f4e59", r=.035, edge="#203746"):
    ax.add_patch(Circle(xy, r, facecolor=color, edgecolor=edge, lw=.8, zorder=5))


def bond(ax, a, b, color=LINE, lw=2.0, ls="-"):
    ax.plot([a[0], b[0]], [a[1], b[1]], color=color, lw=lw, ls=ls,
            solid_capstyle="round", zorder=2)


def arrow(ax, a, b, color=TEAL, lw=1.8, ms=11, style="-|>"):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle=style, mutation_scale=ms,
                                 color=color, lw=lw, zorder=7))


def panel(ax, label, title, subtitle):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.add_patch(FancyBboxPatch((.015, .015), .97, .97,
                                boxstyle="round,pad=.012,rounding_size=.035",
                                facecolor="white", edgecolor="#d9e1e6", lw=1.0))
    ax.text(.055, .91, label, color=TEAL, fontsize=12, weight="bold", va="center")
    ax.text(.13, .91, title, color=NAVY, fontsize=10.5, weight="bold", va="center")
    ax.text(.055, .825, subtitle, color=TEXT, fontsize=7.4, va="center")


def draw(ax, label, title, subtitle, fn):
    panel(ax, label, title, subtitle)
    fn(ax)


def draw_stretch_angle(ax):
    c=np.array([.5,.42]); pts=[np.array([.23,.35]),np.array([.37,.72]),np.array([.70,.70]),np.array([.78,.30])]
    for p in pts: bond(ax,c,p); atom(ax,p,"#63727e",.035)
    atom(ax,c,AMBER,.045)
    arrow(ax, pts[0], pts[0]+np.array([-.11,-.05]), AMBER)
    arrow(ax, pts[2], pts[2]+np.array([.10,.04]), AMBER)
    ax.add_patch(Arc((.5,.42),.32,.27,theta1=35,theta2=145,color=TEAL,lw=2.5))
    ax.text(.5,.15,"stretch + valence angle",ha="center",color=NAVY,fontsize=8)


def draw_salc(ax):
    c=np.array([.5,.48]); angles=np.linspace(0,2*np.pi,6,endpoint=False)+.25
    pts=[c+.20*np.array([np.cos(t),np.sin(t)]) for t in angles]
    for p in pts: bond(ax,c,p,color="#b6c2c9"); atom(ax,p,"#63727e",.033)
    atom(ax,c,AMBER,.043)
    for i,p in enumerate(pts):
        arrow(ax,p,p+.08*(p-c)/np.linalg.norm(p-c), TEAL if i%2==0 else VIOLET,lw=1.5,ms=9)
    ax.text(.5,.12,r"local SALC: $A_1+E$",ha="center",color=NAVY,fontsize=8)


def draw_dihedral(ax):
    pts=[np.array([.20,.35]),np.array([.38,.56]),np.array([.62,.44]),np.array([.80,.65])]
    for a,b in zip(pts[:-1],pts[1:]): bond(ax,a,b)
    for p in pts: atom(ax,p,"#596a76",.038)
    curve=MplPath([pts[1],(.42,.30),(.58,.30),pts[2]],
                  [MplPath.MOVETO,MplPath.CURVE4,MplPath.CURVE4,MplPath.CURVE4])
    ax.add_patch(PathPatch(curve,facecolor="none",edgecolor=VIOLET,lw=2.8,zorder=6))
    ax.text(.5,.18,"one selected exocyclic dihedral",ha="center",color=NAVY,fontsize=8)


def draw_ring(ax):
    ang=np.linspace(0,2*np.pi,7)[:-1]+np.pi/6
    pts=[np.array([.5,.46])+.25*np.array([np.cos(t),np.sin(t)]) for t in ang]
    for a,b in zip(pts,pts[1:]+pts[:1]): bond(ax,a,b)
    for p in pts: atom(ax,p,"#657782",.034)
    for p in pts[:3]: arrow(ax,p,p+.08*(p-np.array([.5,.46])),AMBER,lw=1.4,ms=8)
    ax.add_patch(Arc((.5,.46),.34,.19,theta1=15,theta2=165,color=VIOLET,lw=2.5))
    ax.text(.5,.14,"endocyclic angles • puckering • butterfly",ha="center",color=NAVY,fontsize=7.5)


def draw_oop(ax):
    c=np.array([.5,.43]); pts=[np.array([.25,.36]),np.array([.75,.36]),np.array([.5,.72])]
    for p in pts: bond(ax,c,p); atom(ax,p,"#657782",.035)
    atom(ax,c,AMBER,.043)
    for p in pts: arrow(p if False else ax, p, p) if False else None
    arrow(ax,np.array([.5,.43]),np.array([.5,.18]),RED,lw=2.0)
    ax.text(.5,.12,"signed local SOOP function",ha="center",color=NAVY,fontsize=7.5)


def draw_fragments(ax):
    left=np.array([.29,.49]); right=np.array([.71,.49])
    for center in [left,right]:
        tri=[center+.11*np.array([np.cos(t),np.sin(t)]) for t in (0,2.1,4.2)]
        for a,b in zip(tri,tri[1:]+tri[:1]): bond(ax,a,b)
        for p in tri: atom(ax,p,"#63727e",.028)
    bond(ax,left,right,TEAL,2.0,"--")
    arrow(ax,left+np.array([-.12,-.16]),left+np.array([.03,-.16]),BLUE)
    arrow(ax,right+np.array([.12,.16]),right+np.array([-.03,.16]),VIOLET)
    ax.text(.5,.13,"open graph: pseudobonds  /  closed graph: quaternions",
            ha="center",color=NAVY,fontsize=7.2)


def draw_center(ax):
    ring=[np.array([.33,.38])+ .17*np.array([np.cos(t),np.sin(t)]) for t in np.linspace(0,2*np.pi,7)[:-1]]
    for a,b in zip(ring,ring[1:]+ring[:1]): bond(ax,a,b)
    for p in ring: atom(ax,p,"#657782",.029)
    cen=np.mean(ring,axis=0); ax.scatter(*cen,marker="D",s=80,c=TEAL,edgecolors="white",lw=.8,zorder=6)
    m=np.array([.72,.48]); atom(ax,m,"#9b86b8",.055)
    bond(ax,cen,m,TEAL,2.1,"--")
    arrow(ax,m+np.array([0,.10]),m+np.array([.15,.10]),VIOLET)
    ax.text(.5,.13,r"atom--center and $\eta^n$",ha="center",color=NAVY,fontsize=8)


def draw_bridge(ax):
    a=np.array([.27,.45]); b=np.array([.73,.45]); h=np.array([.5,.58])
    bond(ax,a,h,TEAL,2.2,"--"); bond(ax,h,b,TEAL,2.2,"--")
    atom(ax,a,"#657782",.05); atom(ax,b,"#657782",.05); atom(ax,h,"#e9eef0",.035,"#9aa7ae")
    # heteroatom closure below
    pts=[np.array([.37,.25]),np.array([.63,.25]),np.array([.72,.14]),np.array([.28,.14])]
    for x,y in zip(pts,pts[1:]+pts[:1]): bond(ax,x,y,AMBER,1.8)
    for p in pts: atom(ax,p,"#d58a24",.025)
    ax.text(.5,.72,"electron-deficient bridge",ha="center",color=NAVY,fontsize=8)
    ax.text(.5,.06,"heteroatom ring closures",ha="center",color=NAVY,fontsize=7.5)


def draw_inverse(ax):
    c=np.array([.5,.49]); pts=[c+.21*np.array([np.cos(t),np.sin(t)]) for t in np.linspace(0,2*np.pi,9)[:-1]]
    for i,p in enumerate(pts):
        atom(ax,p,"#88aabd",.035); bond(ax,c,p,"#9aabb5",1.0)
        for q in pts[i+1:]: bond(ax,p,q,"#c6d0d5",.65)
    atom(ax,c,"#477d99",.045)
    ax.text(.5,.13,r"all-pairs $1/R_{ij}$ space",ha="center",color=NAVY,fontsize=8)


def main():
    OUT.mkdir(exist_ok=True)
    fig, axes = plt.subplots(3,3,figsize=(10.2,8.8),facecolor=PALE)
    fig.subplots_adjust(left=.035,right=.965,top=.88,bottom=.045,wspace=.035,hspace=.085)
    fig.suptitle("Coordinate families are chosen from the molecular structure",
                 color=NAVY,fontsize=17,weight="bold",x=.05,ha="left")
    fig.text(.052,.915,"Elementary functions are combined within chemically distinct families before symmetry adaptation.",
             color=TEXT,fontsize=9.5,ha="left")
    entries=[
        ("a","Stretches + angles","bonded valence motions",draw_stretch_angle),
        ("b","SALCs and coordination shells",r"$XY_2/XY_3$, CN 4--14",draw_salc),
        ("c","Exocyclic dihedrals",r"Prelog, then $Z_{\rm eff}$",draw_dihedral),
        ("d","Ring coordinates","angles, torsions, puckering",draw_ring),
        ("e","SOOP functions","endo- and exocyclic motions",draw_oop),
        ("f","Interfragment blocks","open: pseudobonds / closed: quaternions",draw_fragments),
        ("g","Centers and hapticity",r"center--atom, $\eta^n$",draw_center),
        ("h","Bridges and closures","H/halogen and heteroatom cycles",draw_bridge),
        ("i","Topology-independent fallback",r"complete $1/R_{ij}$ pool",draw_inverse),
    ]
    for ax,(lab,title,sub,fn) in zip(axes.flat,entries): draw(ax,lab,title,sub,fn)
    fig.savefig(OUT/"coordinate_family_atlas.pdf",bbox_inches="tight",facecolor=PALE)
    fig.savefig(OUT/"coordinate_family_atlas.png",dpi=360,bbox_inches="tight",facecolor=PALE)


if __name__ == "__main__":
    main()
