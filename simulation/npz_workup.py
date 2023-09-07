import numpy as np
import plotly.graph_objects as go

# Load the data
data = np.load('output_temperatures.npz')
all_T = data['arr_0']
all_T[all_T < 0] = 0  # Set negative values to zero

Nx, Ny, Nz = all_T.shape
x = np.linspace(0, 10, Nx)
y = np.linspace(0, 10, Ny)
z = np.linspace(0, 10, Nz)

# Define the function to get a slice
def get_the_slice(x, y, z, surfacecolor, showscale=False):
    return go.Surface(
        x=x, y=y, z=z,
        surfacecolor=surfacecolor,
        coloraxis='coloraxis',
        showscale=showscale
    )

# Create the initial frame
X, Y = np.meshgrid(x, y)
Z = z[Nz // 2] * np.ones_like(X)
initial_slice = get_the_slice(X, Y, Z, all_T[:, :, Nz // 2], showscale=True)

# Create the animation frames
frames = [go.Frame(
    data=[get_the_slice(X, Y, z[k] * np.ones_like(X), all_T[:, :, k], showscale=(k == 0))],
    name=str(k)
) for k in range(Nz)]

# Define the steps for the animation slider
steps = []
for i in range(Nz):
    step = {
        "args": [
            [str(i)],
            {"frame": {"duration": 300, "redraw": True},
             "mode": "immediate",
             "transition": {"duration": 300}
            }
        ],
        "label": str(i),
        "method": "animate",
    }
    steps.append(step)

# Define the slider
sliders = [{
    "active": 0,
    "yanchor": "top",
    "xanchor": "left",
    "currentvalue": {
        "font": {"size": 20},
        "prefix": "z-slice: ",
        "visible": True,
        "xanchor": "right"
    },
    "transition": {"duration": 300, "easing": "cubic-in-out"},
    "pad": {"b": 10, "t": 50},
    "len": 0.9,
    "x": 0.1,
    "y": 0,
    "steps": steps
}]

# Create the figure
fig = go.Figure(data=[initial_slice], frames=frames)

# Update the layout for the animation
fig.update_layout(
    title="Animation of Temperature Distribution in Z-direction",
    updatemenus=[{
        "buttons": [
            {
                "args": [None, {"frame": {"duration": 500, "redraw": True},
                                "fromcurrent": True, "transition": {"duration": 300, "easing": "quadratic-in-out"}}],
                "label": "Play",
                "method": "animate"
            },
            {
                "args": [[None], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate",
                                  "transition": {"duration": 0}}],
                "label": "Pause",
                "method": "animate"
            }
        ],
        "direction": "left",
        "pad": {"r": 10, "t": 87},
        "showactive": False,
        "type": "buttons",
        "x": 0.1,
        "xanchor": "right",
        "y": 0,
        "yanchor": "top"
    }],
    sliders=sliders,
    scene=dict(zaxis=dict(range=[z[0], z[-1]])),
    coloraxis=dict(colorscale='hot', colorbar_thickness=25, colorbar_len=0.75)
)

fig.show()
