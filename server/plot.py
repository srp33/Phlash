import matplotlib
matplotlib.use("Agg")
from matplotlib.figure import Figure
# from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import pandas as pd


def make_plots(gdata_file, og_start, og_stop, new_starts):
    # matplotlib.rcParams['font.size'] = 8
    gdata_df = pd.read_csv(gdata_file, sep='\t', skiprows=16)
    gdata_df.columns = ['Base', '1', '2', '3', '4', '5', '6']
    gdata_df = gdata_df[gdata_df.Base.isin(range(og_start-100, og_stop+100))]

    fig = Figure()
    axs = fig.subplots(3, sharex=True, sharey=True)
    fig.tight_layout()
    fig.subplots_adjust(left=0.11, bottom=0.13, right=0.96, top=0.87, wspace=0.2, hspace=0.24)
    fig.suptitle('ORF Probabilities')
    for index in range(0, 3):
        frame = str(index + 1)
        axs[index].plot('Base', frame, data=gdata_df)
        axs[index].set_ylim([0, 1])
        axs[index].set_yticks([0, 0.5, 1])
        axs[index].axvline(x=og_start, color='r', linewidth=1)
        axs[index].axvline(x=og_stop, color='r', linewidth=1)
        for new_start in new_starts:
            axs[index].axvline(x=new_start, color='g', linewidth=0.5)
        if index == 0:
            position = og_start - 200
            for new_start in new_starts:
                axs[index].annotate(new_start, 
                    xy=(new_start, 1), 
                    xytext=(position, 1.2), 
                    arrowprops = dict(arrowstyle='->', color='black'))
                position += 150
            axs[index].annotate(og_start, 
                xy=(og_start, 1), 
                xytext=(og_start+100, 1.2), 
                arrowprops = dict(arrowstyle='->', color='black'))
            axs[index].annotate(og_stop, 
                xy=(og_stop, 1), 
                xytext=(og_stop, 1.2), 
                arrowprops = dict(arrowstyle='->', color='black'))
        elif index == 1:
            axs[index].set_ylabel('Direct Sequence')
        else:
            axs[index].set_xlabel('Nucleotide Position')

    # fig.set_size_inches(4, 3)
    fig.savefig(f"../client/public/{og_start}_{og_stop}_direct.png")
    fig.clf()

    axs = fig.subplots(3, sharex=True, sharey=True)
    fig.tight_layout()
    fig.subplots_adjust(left=0.11, bottom=0.13, right=0.96, top=0.91, wspace=0.2, hspace=0.24)
    for index in range(3, 6):
        frame = str(index + 1)
        index = index - 3
        axs[index].plot('Base', frame, data=gdata_df)
        axs[index].set_ylim([0, 1])
        axs[index].set_yticks([0, 0.5, 1])
        axs[index].axvline(x=og_start, color='r', linewidth=1)
        axs[index].axvline(x=og_stop, color='r', linewidth=1)
        for new_start in new_starts:
            axs[index].axvline(x=new_start, color='g', linewidth=0.5)
        if index == 0:
            position = og_start - 200
            for new_start in new_starts:
                axs[index].annotate(new_start, 
                    xy=(new_start, 1), 
                    xytext=(position, 1.2), 
                    arrowprops = dict(arrowstyle='->', color='black'))
                position += 150
            axs[index].annotate(og_start, 
                xy=(og_start, 1), 
                xytext=(og_start+100, 1.2), 
                arrowprops = dict(arrowstyle='->', color='black'))
            axs[index].annotate(og_stop, 
                xy=(og_stop, 1), 
                xytext=(og_stop, 1.2), 
                arrowprops = dict(arrowstyle='->', color='black'))
        elif index == 1:
            axs[index].set_ylabel('Complementary Sequence')
        else:
            axs[index].set_xlabel('Nucleotide Position')

    # fig.set_size_inches(4, 3)
    fig.savefig(f"../client/public/{og_start}_{og_stop}_complementary.png")


# make_plots("uploads/fern.fasta.gdata", 345, 2069, [333, 339])
# make_plot_complementary("uploads/fern.fasta.gdata", 345, 2069, [333, 339])


# def make_plot_complementary(gdata_file, og_start, og_stop, new_starts):
#     gdata_df = pd.read_csv(gdata_file, sep='\t', skiprows=16)
#     gdata_df.columns = ['Base', '1', '2', '3', '4', '5', '6']
#     gdata_df = gdata_df[gdata_df.Base.isin(range(og_start-100, og_stop+100))]

#     fig = Figure()
#     axs = fig.subplots(3, sharex=True, sharey=True)
#     fig.tight_layout()
#     fig.subplots_adjust(left=0.11, bottom=0.13, right=0.97, top=0.87, wspace=0.2, hspace=0.24)
#     for index in range(3, 6):
#         frame = str(index + 1)
#         index = index - 3
#         axs[index].plot('Base', frame, data=gdata_df)
#         axs[index].set_ylim([0, 1])
#         axs[index].set_yticks([0, 0.5, 1])
#         axs[index].axvline(x=og_start, color='r', linewidth=1)
#         axs[index].axvline(x=og_stop, color='r', linewidth=1)
#         for new_start in new_starts:
#             axs[index].axvline(x=new_start, color='g', linewidth=0.5)
#         if index == 0:
#             position = og_start - 200
#             for new_start in new_starts:
#                 axs[index].annotate(new_start, 
#                     xy=(new_start, 1), 
#                     xytext=(position, 1.2), 
#                     arrowprops = dict(arrowstyle='->', color='black'))
#                 position += 150
#             axs[index].annotate(og_start, 
#                 xy=(og_start, 1), 
#                 xytext=(og_start+100, 1.2), 
#                 arrowprops = dict(arrowstyle='->', color='black'))
#             axs[index].annotate(og_stop, 
#                 xy=(og_stop, 1), 
#                 xytext=(og_stop, 1.2), 
#                 arrowprops = dict(arrowstyle='->', color='black'))
#         elif index == 1:
#             axs[index].set_ylabel('Complementary Sequence')
#         else:
#             axs[index].set_xlabel('Nucleotide Position')

#     fig.savefig(f"../client/public/{og_start}_{og_stop}_complementary.png")

# gdata_file = "uploads/fern.fasta.gdata"
# gdata_df = pd.read_csv(gdata_file, sep='\t', skiprows=16)
# gdata_df.columns = ['Base', '1', '2', '3', '4', '5', '6']

