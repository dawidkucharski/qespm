"""Regenerate manuscript/figT1_traceability.pdf.

Measurement traceability chain of QESPM:
  SI metre  ->  stage/scanner calibration  ->  ion-surface separation h
  SI second ->  secular frequency (omega_x, Delta omega_x)
  sample surface  --sensed by the ion-->  Delta omega_x
  h  --calibrates-->  calibrated forward model (ITF)
  ITF  --validates-->  regularised reconstruction
  Delta omega_x  --inverted by-->  reconstruction  --estimates-->  z_s(x,y)
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

EDGE = '#1f4e79'
FILL = '#dce9f5'
TEXT = '#0b1f33'

fig, ax = plt.subplots(figsize=(12.6, 7.8))
ax.set_xlim(0, 12.6)
ax.set_ylim(-2.3, 8.1)
ax.axis('off')

BOX = dict(boxstyle='round,pad=0.5,rounding_size=0.4',
           facecolor=FILL, edgecolor=EDGE, linewidth=1.6)


def box(x, y, w, h, text, fs=13):
    ax.add_patch(FancyBboxPatch((x, y), w, h, **BOX))
    ax.text(x + w / 2, y + h / 2, text, ha='center', va='center',
            fontsize=fs, color=TEXT)


def arrow(p1, p2, label=None, lpos=None):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle='-|>',
                                 mutation_scale=22, color=EDGE,
                                 linewidth=1.8, connectionstyle='arc3,rad=0'))
    if label:
        ax.text(*lpos, label, ha='center', va='center', fontsize=11.5,
                color=EDGE, style='italic',
                bbox=dict(facecolor='white', edgecolor='none', pad=1.6))


# --- boxes (x, y, w, h) ---
A1 = (0.4, 5.6)   # SI metre
A2 = (0.4, 3.1)   # stage / scanner
A3 = (0.4, 0.6)   # h
B1 = (8.6, 5.6)   # SI second
B2 = (8.6, 3.1)   # secular frequency
C1 = (4.4, 3.1)   # sample surface
F1 = (4.4, 0.6)   # forward model (ITF)
D1 = (4.4, -1.9)  # reconstruction
D2 = (8.6, -1.9)  # measurand
W, H = 3.4, 1.8

box(*A1, W, H, 'SI metre\n(calibrated grating)')
box(*A2, W, H, 'stage / scanner\ncalibration')
box(*A3, W, H, 'ion\u2013surface\nseparation $h$')
box(*B1, W, H, 'SI second\n(frequency reference)')
box(*B2, W, H, 'secular frequency\n$\\omega_x,\\ \\Delta\\omega_x$')
box(*C1, W, H, 'sample surface\n(topography, charge)')
box(*F1, W, H, 'calibrated forward\nmodel (ITF)')
box(*D1, W, H, 'regularised\nreconstruction')
box(*D2, W, H, '$z_s(x,y)$\n(measurand)')

# --- arrows ---
arrow((2.1, 5.6), (2.1, 4.9), 'calibrates', (2.85, 5.25))
arrow((2.1, 3.1), (2.1, 2.4), 'calibrates', (2.85, 2.75))
arrow((3.8, 1.5), (4.4, 1.5), 'calibrates', (4.1, 1.82))
arrow((10.3, 5.6), (10.3, 4.9), 'calibrates', (11.05, 5.25))
arrow((7.8, 4.0), (8.6, 4.0), 'senses', (8.2, 4.35))
arrow((6.1, 0.6), (6.1, -0.1), 'validates', (6.95, 0.25))
arrow((9.2, 3.1), (7.6, -0.1), 'inverts', (8.6, 1.5))
arrow((7.8, -1.0), (8.6, -1.0), 'estimates', (8.2, -0.62))

fig.savefig('manuscript/figT1_traceability.pdf', bbox_inches='tight',
            pad_inches=0.12)
print('saved manuscript/figT1_traceability.pdf')
