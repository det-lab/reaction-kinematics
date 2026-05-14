import pandas as pd

from reaction_kinematics.masstable import MTAB

rows = []

for isotope, mass in MTAB.items():
    rows.append({
        "Isotope": isotope,
        "Mass (MeV/c²)": mass,
    })

df = pd.DataFrame(rows).sort_values("Isotope")

html = df.to_html(
    index=False,
    table_id="mass-table",
    classes=["display", "mass-table"]
)

with open("docs/includes/mass_table.html", "w") as f:
    f.write(html)
    