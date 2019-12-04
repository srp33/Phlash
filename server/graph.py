import altair as alt
import pandas as pd
from altair import datum
from statistics import mean
# from IPython.display import display


def make_graph_direct(df, og_start, og_stop, new_starts):
    """
    Creates graphs containing direct sequences (frames 1, 2, 3).
    :return: altair chart as a dictionary
    """
    # Create a selection that chooses the nearest point & selects based on x-value
    nearest = alt.selection_single(nearest=True, on='mouseover', empty='none')
    
    # frame 1
    one = alt.Chart(df).mark_line().encode(
        x=alt.X('Base', axis=alt.Axis(title='')),
        y=alt.Y('1', axis=alt.Axis(title=''), scale=alt.Scale(domain=(0, 1)))
    )

    # frame 2
    two = alt.Chart(df).mark_line().encode(
        x=alt.X('Base', axis=alt.Axis(title='')),
        y=alt.Y('2', axis=alt.Axis(title='Direct Sequences'), scale=alt.Scale(domain=(0, 1)))
    )

    # frame 3
    three = alt.Chart(df).mark_line().encode(
        x=alt.X('Base', axis=alt.Axis(title='Nucleotide Position')),
        y=alt.Y('3', axis=alt.Axis(title=''), scale=alt.Scale(domain=(0, 1)))
    )

    # Transparent selectors across the chart. This is what tells us the y-value of the cursor.
    selectors = alt.Chart().mark_point().encode(
        x='Base:Q',
        opacity=alt.value(0),
    ).add_selection(
        nearest
    )

    # Draw points on the line, and highlight based on selection
    points = one.mark_point().encode(
        opacity=alt.condition(nearest, alt.value(1), alt.value(0))
    )

    # Draw text labels near the points, and highlight based on selection
    text = one.mark_text(align='left', dx=5, dy=-10).encode(
        text=alt.condition(nearest, 'Base:Q', alt.value(' '))
    )

    # Draw a rule at the location of the selection
    rules = alt.Chart().mark_rule(color='grey').encode(
        x='Base:Q',
    ).transform_filter(
        nearest
    )

    # vline at original start of DNA Master gene call
    start = alt.Chart(df).mark_rule(color='red').encode(
        x='a:Q',
        size=alt.value(0.5)
    ).transform_calculate(
        a=str(og_start)
    )

    # vline at original stop of DNA Master gene call
    stop = alt.Chart(df).mark_rule(color='red').encode(
        x='a:Q',
        size=alt.value(0.5)
    ).transform_calculate(
        a=str(og_stop)
    )

    # Put the layers into a chart and bind the data
    one_layered = alt.layer(
        one, selectors, points, rules, text, start, stop, data=df, width=600, height=150
    )

    two_layered = alt.layer(
        two, selectors, points, rules, text, start, stop, data=df, width=600, height=150
    )

    three_layered = alt.layer(
        three, selectors, points, rules, text, start, stop, data=df, width=600, height=150
    )

    chart = alt.vconcat(
        one_layered, two_layered, three_layered
    ).configure_axis(
        labelFontSize=12,
        titleFontSize=13
    )

    return chart.to_dict()


def make_graph_complementary(df, og_start, og_stop, new_starts):
    """
    Creates graphs containing complementary sequences (frames 4, 5, 6).
    :return: altair chart as a dictionary
    """
    # Create a selection that chooses the nearest point & selects based on x-value
    nearest = alt.selection_single(nearest=True, on='mouseover', empty='none')
    
    # frame 4
    four = alt.Chart(df).mark_line().encode(
        x=alt.X('Base', axis=alt.Axis(title='')),
        y=alt.Y('4', axis=alt.Axis(title=''), scale=alt.Scale(domain=(0, 1)))
    )

    # frame 5
    five = alt.Chart(df).mark_line().encode(
        x=alt.X('Base', axis=alt.Axis(title='')),
        y=alt.Y('5', axis=alt.Axis(title='Complementary Sequences'), scale=alt.Scale(domain=(0, 1)))
    )

    # frame 6
    six = alt.Chart(df).mark_line().encode(
        x=alt.X('Base', axis=alt.Axis(title='Nucleotide Position')),
        y=alt.Y('6', axis=alt.Axis(title=''), scale=alt.Scale(domain=(0, 1)))
    )

    # Transparent selectors across the chart. This is what tells us the y-value of the cursor.
    selectors = alt.Chart().mark_point().encode(
        x='Base:Q',
        opacity=alt.value(0),
    ).add_selection(
        nearest
    )

    # Draw points on the line, and highlight based on selection
    points = four.mark_point().encode(
        opacity=alt.condition(nearest, alt.value(1), alt.value(0))
    )

    # Draw text labels near the points, and highlight based on selection
    text = four.mark_text(align='left', dx=5, dy=-10).encode(
        text=alt.condition(nearest, 'Base:Q', alt.value(' '))
    )

    # Draw a rule at the location of the selection
    rules = alt.Chart().mark_rule(color='gray').encode(
        x='Base:Q',
    ).transform_filter(
        nearest
    )

    # vline at original start of DNA Master gene call
    start = alt.Chart(df).mark_rule(color='red').encode(
        x='a:Q',
        size=alt.value(0.5)
    ).transform_calculate(
        a=str(og_start)
    )

    # vline at original stop of DNA Master gene call
    stop = alt.Chart(df).mark_rule(color='red').encode(
        x='a:Q',
        size=alt.value(0.5)
    ).transform_calculate(
        a=str(og_stop)
    )

    # Put the layers into a chart and bind the data
    one_layered = alt.layer(
        four, selectors, points, rules, text, start, stop, data=df, width=600, height=150
    )

    two_layered = alt.layer(
        five, selectors, points, rules, text, start, stop, data=df, width=600, height=150
    )

    three_layered = alt.layer(
        six, selectors, points, rules, text, start, stop, data=df, width=600, height=150
    )

    chart = alt.vconcat(
        one_layered, two_layered, three_layered
    ).configure_axis(
        labelFontSize=12,
        titleFontSize=13
    )

    return chart.to_dict()


#--------------------------------------------------------------------------------------
    # .properties(
    #     title='ORF Coding Potentials'
    # ).configure_title(
    #     fontSize=15,
    #     anchor='middle'

# df = pd.read_csv("uploads/fern.fasta.gdata", sep='\t', skiprows=16)
# df.columns = ['Base', '1', '2', '3', '4', '5', '6']

# start = 3303
# stop = 4025
# new_starts = [333, 339]

# make_graph_direct(df, start, stop, new_starts)

# # helper functions
# def get_avg_probs(df, num):
#     probs = []
#     df_sort = df.loc[(df['Base']-num).abs().argsort()[:2]]
#     print(df_sort)
#     for frame in range(1,7):
#         frame = str(frame)
#         prob = mean(df_sort[frame])
#         probs.append(prob)
#     return probs

# def check_base(df, num):
#     # for num in nums:
#     if num not in df.Base.values:
#         probs = get_avg_probs(df, num)
#         print(probs)
#         df = df.append({'Base': num,
#                 '1': probs[0],
#                 '2': probs[1],
#                 '3': probs[2],
#                 '4': probs[3],
#                 '5': probs[4],
#                 '6': probs[5]},
#                 ignore_index=True)
#     return df


