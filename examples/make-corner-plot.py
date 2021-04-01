
# %%
import corner
import numpy as np

# %%
chains = np.loadtxt("chains-2.txt")

# %%

#%%
figure = corner.corner(
    chains,
    labels=[
        "f",
        "a - au",
        "τ - days",
        "inc - °",
        "$\Omega$ - °",
        "ecc",
        "$\omega$ - °",
    ],
    quantiles=[0.16, 0.5, 0.84],
    show_titles=True, title_kwargs={"fontsize":12},
)
# %%
figure.savefig("temp-corner-3.svg", dpi=200)
# %%
