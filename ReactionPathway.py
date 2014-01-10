'''
Created on Jul 7, 2013

@author: Daniel
'''

class ReactionPathway(object):

    def __init__(self, mgi_id, reactome_id, reaction_name):
        self.mgi_id = mgi_id
        self.reactome_id = reactome_id
        self.reaction_name = reaction_name
