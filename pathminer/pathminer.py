# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import logging

from collections import defaultdict
"""
segment mining for BioCyc segments

Takes a series of metabolites, genes or proteins and returns a list of segments and/or
reactions that the most altered in definable ways.

Supply data as a list of tuples e.g. [(BIOCYC-OBJECT, VALUE)]

Supply multiple datasets as list of lists.

Martin Fitzpatrick
May 2014

"""

TARGET_PATHWAYS = 0
TARGET_REACTIONS = 1
TARGET_COMPARTMENTS = 2


def mining(data, target=TARGET_PATHWAYS, algorithm='c', no_of_results=5, shared=True, relative=False):
    
    # Iterate all the compounds in the current analysis
    # Assign score to each of the compound's segments
    # Sum up, crop and return a list of segment_ids to display
    # Pass this in as the list to view
    # + requested segments, - excluded segments

    segment_scores = defaultdict(int)
    logging.info("Mining %s using '%s'" % (target, algorithm))

    if type(data) != list:
        data = [data]

    # Data provided as a list of tuples (ENTITY, data)
    for n, entity_data in enumerate( data ):
        if entity_data is None:
            continue
        entity, score = entity_data

        if entity is None or score is None:
            continue

        # Get the entity's segments
        try:
            if target == TARGET_PATHWAYS:
                segments = entity.pathways
        
            elif target == TARGET_REACTIONS:
                segments = entity.reactions
            
            elif target == TARGET_COMPARTMENTS:
                segments = entity.compartments
        except:
            continue

        if segments == []:
            continue

        if shared:
            # Share the change score between the associated segments
            # this prevents compounds having undue influence
            score = float(score) / len(segments)
        

        for p in segments:
            mining_val = {
                'c': abs(score),
                'u': max(0, score),
                'd': abs(min(0, score)),
                'm': 1.0,
                't': score,
                }
            segment_scores[p] += mining_val[algorithm]


    # If we're using tendency scaling; abs the scores here
    if algorithm == 't':
        for p, v  in list(segment_scores.items()):
            segment_scores[p] = abs(v)


    # If we're pruning, then remove any segments not in keep_segments
    if relative:
        logging.info("Scaling scores to target sizes...")
        for p, v in list(segment_scores.items()):
            segment_scores[p] = float(v) / len(p.reactions)

    # Check if we've got any data
    assert segment_scores, "Not enough data to do anything useful. Add more data, or change the mining type."

    # Now take the accumulated scores; and create the output
    segment_scorest = list(segment_scores.items())  # Switch it to a dict so we can sort
    segment_scorest = [(p, v) for p, v in segment_scorest if v > 0]  # Remove any scores of 0
    segment_scorest.sort(key=lambda tup: tup[1], reverse=True)  # Sort by scores (either system)

    # Get top N defined by no_of_results parameter
    keep_segments = segment_scorest[0:no_of_results]
    remaining_segments = segment_scorest[no_of_results + 1:no_of_results + 100]

    logging.info("Mining recommended %d out of %d" % (len(keep_segments), len(segment_scores)))

    for n, p in enumerate(keep_segments):
        logging.info("- %d. %s [%.2f]" % (n + 1, p[0].name, p[1]))


    return keep_segments
    