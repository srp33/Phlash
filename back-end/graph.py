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


# def make_graph_complementary(df, og_start, og_stop, new_starts):
#     """
#     Creates graphs containing complementary sequences (frames 4, 5, 6).
#     :return: altair chart as a dictionary
#     """
#     # Create a selection that chooses the nearest point & selects based on x-value
#     nearest = alt.selection_single(nearest=True, on='mouseover', empty='none')
    
#     # frame 4
#     four = alt.Chart(df).mark_line().encode(
#         x=alt.X('Base', axis=alt.Axis(title='')),
#         y=alt.Y('4', axis=alt.Axis(title=''), scale=alt.Scale(domain=(0, 1)))
#     )

#     # frame 5
#     five = alt.Chart(df).mark_line().encode(
#         x=alt.X('Base', axis=alt.Axis(title='')),
#         y=alt.Y('5', axis=alt.Axis(title='Complementary Sequences'), scale=alt.Scale(domain=(0, 1)))
#     )

#     # frame 6
#     six = alt.Chart(df).mark_line().encode(
#         x=alt.X('Base', axis=alt.Axis(title='Nucleotide Position')),
#         y=alt.Y('6', axis=alt.Axis(title=''), scale=alt.Scale(domain=(0, 1)))
#     )

#     # Transparent selectors across the chart. This is what tells us the y-value of the cursor.
#     selectors = alt.Chart().mark_point().encode(
#         x='Base:Q',
#         opacity=alt.value(0),
#     ).add_selection(
#         nearest
#     )

#     # Draw points on the line, and highlight based on selection
#     points = four.mark_point().encode(
#         opacity=alt.condition(nearest, alt.value(1), alt.value(0))
#     )

#     # Draw text labels near the points, and highlight based on selection
#     text = four.mark_text(align='left', dx=5, dy=-10).encode(
#         text=alt.condition(nearest, 'Base:Q', alt.value(' '))
#     )

#     # Draw a rule at the location of the selection
#     rules = alt.Chart().mark_rule(color='gray').encode(
#         x='Base:Q',
#     ).transform_filter(
#         nearest
#     )

#     # vline at original start of DNA Master gene call
#     start = alt.Chart(df).mark_rule(color='red').encode(
#         x='a:Q',
#         size=alt.value(0.5)
#     ).transform_calculate(
#         a=str(og_start)
#     )

#     # vline at original stop of DNA Master gene call
#     stop = alt.Chart(df).mark_rule(color='red').encode(
#         x='a:Q',
#         size=alt.value(0.5)
#     ).transform_calculate(
#         a=str(og_stop)
#     )

#     # Put the layers into a chart and bind the data
#     one_layered = alt.layer(
#         four, selectors, points, rules, text, start, stop, data=df, width=600, height=150
#     )

#     two_layered = alt.layer(
#         five, selectors, points, rules, text, start, stop, data=df, width=600, height=150
#     )

#     three_layered = alt.layer(
#         six, selectors, points, rules, text, start, stop, data=df, width=600, height=150
#     )

#     chart = alt.vconcat(
#         one_layered, two_layered, three_layered
#     ).configure_axis(
#         labelFontSize=12,
#         titleFontSize=13
#     )

#     return chart.to_dict()


#--------------------------------------------------------------------------------------
def apex_direct(genemark_gdata_file, start, stop, start_options):
   gdata_df = pd.read_csv(genemark_gdata_file, sep='\t', skiprows=16)
   gdata_df.columns = ['Base', '1', '2', '3', '4', '5', '6']
   gdata_df = gdata_df[gdata_df.Base.isin(range(start-100, stop+100))]
   xaxis = gdata_df["Base"].to_list()
   frame_1 = gdata_df["2"].to_list()

   print(frame_1)


# apex_direct("uploads/fern.fasta.gdata", 345, 2069, [333, 339])

def merge(list1, list2): 
   merged_list = [] 
   for i in range(max((len(list1), len(list2)))): 
      while True: 
         try: 
            tup = [list1[i], list2[i]]
         except IndexError: 
            if len(list1) > len(list2): 
               list2.append('') 
               tup = [list1[i], list2[i]]
            elif len(list1) < len(list2): 
               list1.append('') 
               tup = [list1[i], list2[i]]
            continue
         merged_list.append(tup) 
         break
   return merged_list 
  