#!/usr/bin/env python3
"""Generate the quantitative QESPM sensitivity--standoff landscape."""

import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt

fig, ax1 = plt.subplots(figsize=(8.6, 6.0))

methods = {
    'AFM (contact)':        {'lateral': 1.0,     'vertical': 0.01,   'standoff': 0.3,    'color': '#607D8B', 'marker': 's'},
    'AFM (tapping)':        {'lateral': 5.0,     'vertical': 0.05,   'standoff': 5.0,    'color': '#78909C', 'marker': 's'},
    'STM':                  {'lateral': 0.1,     'vertical': 0.001,  'standoff': 0.5,    'color': '#455A64', 'marker': 'D'},
    'Optical profilometry': {'lateral': 500.0,   'vertical': 1.0,    'standoff': 1e5,    'color': '#FF9800', 'marker': 'o'},
    'SEM':                  {'lateral': 2.0,     'vertical': 1.0,    'standoff': 1e4,    'color': '#9E9E9E', 'marker': '^'},
    r'$\mathbf{QESPM}$':    {'lateral': 20.0,    'vertical': 0.007,  'standoff': 4e4,    'color': '#C62828', 'marker': '*'},
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
        'SEM':                  (6.0, 1.2),
        r'$\mathbf{QESPM}$':    (0.5, 4.0),
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
ax1.set_title('Detection sensitivity--standoff landscape', fontsize=12, fontweight='bold')

# Wide enough to show all techniques and the reconstruction-uncertainty marker
ax1.set_xlim(0.08, 3e5)
ax1.set_ylim(0.0003, 1200)
ax1.grid(True, alpha=0.25, which='both')

# Shade QESPM quadrant (1--100 um standoff, per the operating envelope)
ax1.axvspan(1e3, 1e5, alpha=0.06, color='red')
ax1.axhspan(0.0003, 0.1, alpha=0.06, color='red')
ax1.annotate('QESPM\nniche', xy=(2e4, 0.0008), fontsize=10, color='darkred',
             fontweight='bold', fontstyle='italic')

# Reconstruction uncertainty shown SEPARATELY: it is not part of the
# sensitivity scale (charge-limited, Table: uncertainty budget, h = 40 um)
ax1.plot([4e4, 4e4], [0.007, 319], ls='--', color='#C62828', lw=1.1, alpha=0.6)
ax1.scatter(4e4, 319, marker='o', facecolors='none', edgecolors='#C62828',
            s=170, zorder=10, linewidths=1.7)
ax1.annotate('QESPM reconstruction uncertainty\n(charge-limited, $u_c(z)$ = 319 nm at $h$ = 40 $\mu$m)',
             xy=(4e4, 319), xytext=(1.6e4, 30),
             fontsize=8, color='#C62828',
             arrowprops=dict(arrowstyle='->', color='#C62828', lw=0.9))

plt.tight_layout()
plt.savefig('manuscript/fig4_landscape.pdf', dpi=150, bbox_inches='tight')
plt.close()
print('Saved: manuscript/fig4_landscape.pdf')
