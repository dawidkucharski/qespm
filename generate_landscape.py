#!/usr/bin/env python3
"""
Generate QESPM positioning figure: resolution--sensitivity--standoff landscape.
Clean, non-overlapping labels, all data visible within axes.
"""

import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

# ============================================================
# Panel (a): Resolution--standoff landscape
# ============================================================
methods = {
    'AFM (contact)':        {'lateral': 1.0,     'vertical': 0.01,   'standoff': 0.3,    'color': '#607D8B', 'marker': 's'},
    'AFM (tapping)':        {'lateral': 5.0,     'vertical': 0.05,   'standoff': 5.0,    'color': '#78909C', 'marker': 's'},
    'STM':                  {'lateral': 0.1,     'vertical': 0.001,  'standoff': 0.5,    'color': '#455A64', 'marker': 'D'},
    'Optical profilometry': {'lateral': 500.0,   'vertical': 1.0,    'standoff': 1e5,    'color': '#FF9800', 'marker': 'o'},
    'Confocal microscopy':  {'lateral': 200.0,   'vertical': 5.0,    'standoff': 5e4,    'color': '#FFB74D', 'marker': 'o'},
    'SEM':                  {'lateral': 2.0,     'vertical': 2.0,    'standoff': 1e4,    'color': '#9E9E9E', 'marker': '^'},
    'NV centre (magnetic)': {'lateral': 10.0,    'vertical': None,   'standoff': 50.0,   'color': '#7B1FA2', 'marker': 'P'},
    r'$\mathbf{QESPM}$':    {'lateral': 20.0,    'vertical': 0.01,   'standoff': 1e4,    'color': '#C62828', 'marker': '*'},
}

for name, m in methods.items():
    if m['vertical'] is not None:
        ax1.scatter(m['standoff'], m['vertical'], c=m['color'], s=250 if 'QESPM' in name else 120,
                    marker=m['marker'], zorder=10 if 'QESPM' in name else 5,
                    edgecolors='black' if 'QESPM' in name else 'white', linewidths=1.5 if 'QESPM' in name else 0.5)
    
    # Position labels in data coordinates — manually placed to avoid overlap
    label_x = m['standoff']
    label_y = m['vertical'] if m['vertical'] is not None else 0.01
    
    # Manually adjust each label position
    offsets = {
        'AFM (contact)':        (1.5, 2.0),
        'AFM (tapping)':        (1.8, 1.8),
        'STM':                  (1.5, 2.5),
        'Optical profilometry':  (0.25, 0.5),
        'Confocal microscopy':   (0.3, 0.4),
        'SEM':                  (2.5, 0.5),
        'NV centre (magnetic)': (0.4, 2.5),
        r'$\mathbf{QESPM}$':    (0.2, 3.0),
    }
    mx, my = offsets.get(name, (1.5, 1.5))
    
    if m['vertical'] is not None:
        ax1.annotate(name.replace(r'\mathbf{', '').replace('}', ''),
                     (label_x, label_y),
                     xytext=(label_x * mx, label_y * my),
                     fontsize=8, color=m['color'],
                     fontweight='bold' if 'QESPM' in name else 'normal',
                     arrowprops=dict(arrowstyle='-', color=m['color'], lw=0.5, alpha=0.5) if 'QESPM' in name else None)

ax1.set_xscale('log'); ax1.set_yscale('log')
ax1.set_xlabel('Standoff distance [nm]', fontsize=11)
ax1.set_ylabel('Intrinsic vertical detection sensitivity [nm]', fontsize=11)
ax1.set_title('(a) Detection sensitivity--standoff landscape', fontsize=12, fontweight='bold')

# Wide enough to show all techniques and the reconstruction-uncertainty marker
ax1.set_xlim(0.08, 3e5)
ax1.set_ylim(0.0003, 1200)
ax1.grid(True, alpha=0.25, which='both')

# Shade QESPM quadrant
ax1.axvspan(500, 5e4, alpha=0.06, color='red')
ax1.axhspan(0.0003, 0.1, alpha=0.06, color='red')
ax1.annotate('QESPM\nniche', xy=(3000, 0.002), fontsize=10, color='darkred',
             fontweight='bold', fontstyle='italic')

# Reconstruction uncertainty shown SEPARATELY: it is not part of the
# sensitivity scale (charge-limited, Table: uncertainty budget)
ax1.plot([1e4, 1e4], [0.01, 319], ls='--', color='#C62828', lw=1.1, alpha=0.6)
ax1.scatter(1e4, 319, marker='o', facecolors='none', edgecolors='#C62828',
            s=170, zorder=10, linewidths=1.7)
ax1.annotate('QESPM reconstruction uncertainty\n(charge-limited, $u_c(z)$ = 319 nm)',
             xy=(1e4, 319), xytext=(9e3, 12),
             fontsize=8, color='#C62828',
             arrowprops=dict(arrowstyle='->', color='#C62828', lw=0.9))

# ============================================================
# Panel (b): Radar/spider chart — capability dimensions
# ============================================================
categories = ['Vertical\nsensitivity', 'Lateral\nresolution', 'Standoff\nrange',
              'UHV\ncompatibility', 'Cryogenic\noperation', 'Non-contact', 'Traceable\ncalibration']
N = len(categories)
angles = np.linspace(0, 2*np.pi, N, endpoint=False).tolist()
angles += angles[:1]

# Current state-of-the-art (grey)
sota = [4, 5, 2, 1, 1, 3, 5]
sota += sota[:1]

# QESPM (red)
qespm = [5, 2, 5, 5, 5, 5, 4]
qespm += qespm[:1]

ax2 = plt.subplot(122, polar=True)
ax2.fill(angles, sota, alpha=0.2, color='grey', label='State of the art')
ax2.plot(angles, sota, 'o-', color='grey', lw=2, markersize=6)
ax2.fill(angles, qespm, alpha=0.25, color='red', label='QESPM')
ax2.plot(angles, qespm, 'o-', color='red', lw=2.5, markersize=8)

ax2.set_xticks(angles[:-1])
ax2.set_xticklabels(categories, fontsize=9)
ax2.set_ylim(0, 5.5)
ax2.set_yticks([1, 2, 3, 4, 5])
ax2.set_yticklabels(['1', '2', '3', '4', '5'], fontsize=7)
ax2.set_title('(b) Capability profile (1–5 scale)', fontsize=12, fontweight='bold', pad=20)
ax2.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1), fontsize=9)

plt.tight_layout()
plt.savefig('manuscript/fig4_landscape.pdf', dpi=150, bbox_inches='tight')
plt.close()
print('Saved: manuscript/fig4_landscape.pdf')
