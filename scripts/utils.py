from pathlib import Path

FIGDIR = Path("figures")


def save_fig(ax, filename):
    """Save the figure behind a scanpy plot (Axes, list of Axes, or Figure) to figures/."""
    fig = ax[0].figure if isinstance(ax, list) else getattr(ax, "figure", ax)
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGDIR / filename, dpi=150, bbox_inches="tight")
