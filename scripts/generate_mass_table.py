import pandas as pd

from reaction_kinematics.masstable import MTAB

rows = []

for (A, elem), mass in MTAB.items():
    isotope = f"{A}{elem.capitalize()}"  # -> 3He, 12C, etc.

    rows.append({
        "Isotope": isotope,
        "Mass (MeV/c²)": mass,
    })

df = pd.DataFrame(rows).sort_values("Isotope")

table_md = df.to_markdown(index=False)

# Full markdown page content
markdown = f"""# Mass Table

On this page you can find the mass table this library uses.

{table_md}
"""

with open("../docs/mass_table.md", "w", encoding="utf-8") as f:
    f.write(markdown)