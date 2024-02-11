from ascii_graph import Pyasciigraph
from ascii_graph.colors import *
from ascii_graph.colordata import vcolor
from ascii_graph.colordata import hcolor

test = [('long_label', 423), ('sl', 1234), ('line3', 531),
    ('line4', 200), ('line5', 834)]

# One color per line
print('Color example:')
pattern = [Gre, Yel, Red]
data = vcolor(test, pattern)

graph = Pyasciigraph()
for line in graph.graph('vcolor test', data):
    print(line)

# Multicolor on one line
print('\nMultiColor example:')

# Color lines according to Thresholds
thresholds = {
  51:  Gre, 100: Blu, 350: Yel, 500: Red,
}
data = hcolor(test, thresholds)

# graph with colors, power of 1000, different graph symbol,
# float formatting and a few tweaks
graph = Pyasciigraph(
    line_length=120,
    min_graph_length=50,
    separator_length=4,
    multivalue=False,
    human_readable='si',
    graphsymbol='*',
    float_format='{0:,.2f}',
    force_max_value=2000,
    )

for line in graph.graph(label=None, data=data):
    print(line)