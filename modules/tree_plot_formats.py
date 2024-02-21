class TreeLineFormats(object):
    @classmethod
    def color_by_binary_trait(node):
        if node.trait == 1:
            color = 'red'
        elif node.trait == 0:
            color = 'blue'
        else:
            raise Exception('Unknown trait value:', node.trait)
        return {'color': color}