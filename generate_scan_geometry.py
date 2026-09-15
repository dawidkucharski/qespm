#!/usr/bin/env python3
"""Generate the QESPM scan geometry schematic."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyArrowPatch, FancyBboxPatch


EDGE = '#1f4e79'
FILL = '#dce9f5'
ION = '#b71c1c'
SAMPLE = '#8d6e63'
TEXT = '#0b1f33'


fig, ax = plt.subplots(figsize=(8.2, 4.1))
ax.set_xlim(0, 10)
ax.set_ylim(0, 5.2)
ax.axis('off')

trap = FancyBboxPatch((1.2, 4.05), 7.6, 0.6,
                      boxstyle='round,pad=0.15,rounding_size=0.18',
                      facecolor=FILL, edgecolor=EDGE, linewidth=1.8)
ax.add_patch(trap)
ax.add_patch(Rectangle((3.1, 3.77), 0.9, 0.18, facecolor=EDGE, edgecolor=EDGE))
ax.add_patch(Rectangle((6.0, 3.77), 0.9, 0.18, facecolor=EDGE, edgecolor=EDGE))
ax.text(5.0, 4.37, 'planar trap / electrodes', ha='center', va='center',
        fontsize=11.5, color=TEXT)

ion_xy = (5.0, 2.95)
ax.add_patch(Circle(ion_xy, 0.17, facecolor=ION, edgecolor='black', linewidth=1.0))
ax.text(5.0, 3.46, 'stationary trapped ion', ha='center', va='center',
        fontsize=10.5, color=ION)

ax.add_patch(FancyArrowPatch((5.0, 2.75), (5.0, 1.95), arrowstyle='-|>',
                             mutation_scale=18, linewidth=1.8, color=EDGE))
ax.text(5.23, 2.28, '$\\Delta\\omega_x(x,y)$', ha='left', va='center',
        fontsize=12, color=EDGE)

sample = Rectangle((1.5, 0.9), 7.0, 0.48, facecolor='#d7ccc8', edgecolor=SAMPLE, linewidth=1.4)
ax.add_patch(sample)
ax.plot([1.5, 2.1, 2.7, 3.5, 4.2, 5.0, 5.7, 6.5, 7.2, 8.0, 8.5],
        [1.38, 1.46, 1.31, 1.44, 1.28, 1.41, 1.35, 1.47, 1.30, 1.42, 1.36],
        color=SAMPLE, linewidth=2.0)
ax.text(5.0, 0.62, 'translated sample surface', ha='center', va='center',
        fontsize=12, color=TEXT)

ax.add_patch(FancyArrowPatch((2.3, 0.42), (7.7, 0.42), arrowstyle='<->',
                             mutation_scale=16, linewidth=1.6, color=EDGE))
ax.text(5.0, 0.04, 'lateral scan $x,y$', ha='center', va='bottom',
        fontsize=11, color=EDGE)

ax.plot([4.25, 4.25], [1.42, 2.8], color='black', linestyle='--', linewidth=1.1)
ax.add_patch(FancyArrowPatch((4.25, 2.75), (4.25, 1.48), arrowstyle='<->',
                             mutation_scale=14, linewidth=1.2, color='black'))
ax.text(4.02, 2.23, '$h$', ha='right', va='center', fontsize=13, color='black')

note = FancyBboxPatch((7.0, 1.95), 2.35, 0.78,
                      boxstyle='round,pad=0.18,rounding_size=0.12',
                      facecolor='white', edgecolor=EDGE, linewidth=1.0)
ax.add_patch(note)
ax.text(8.175, 2.34, 'local region\napproximates grounded\nhalf-space',
        ha='center', va='center', fontsize=10, color=TEXT)

fig.savefig('manuscript/figS2_scan.pdf', bbox_inches='tight', pad_inches=0.06)
fig.savefig('figS2_scan.pdf', bbox_inches='tight', pad_inches=0.06)
print('saved manuscript/figS2_scan.pdf and figS2_scan.pdf')