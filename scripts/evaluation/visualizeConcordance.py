import pandas as pd
import altair as alt

df = pd.read_csv(snakemake.input[0])
# Define additional columns
def classification(row):
    if row.illuminaDropout or row.nanoporeDropout:
        return "Dropout"
    elif row.concordance:
        return "Concordance"
    else:
        return "Discordance"


df["classification"] = df.apply(lambda row: classification(row), axis=1)

# Sample Selector
make = pd.DataFrame({"sample": list(df["sample"].unique())})
selection_sm = alt.selection_multi(fields=["sample"])
color = alt.condition(selection_sm, alt.value("green"), alt.value("lightgray"))
# Add Boolean Selectors
selector_fn = alt.binding_select(options=[True, False, None], name="False Negatives")
selection_fn = alt.selection_single(fields=["falseNegative"], bind=selector_fn)
selector_fp = alt.binding_select(options=[True, False, None], name="False Positives")
selection_fp = alt.selection_single(fields=["falsePositive"], bind=selector_fp)
selector_dh = alt.binding_select(options=[True, False, None], name="Detected HETs")
selection_dh = alt.selection_single(fields=["detectedHET"], bind=selector_dh)
selector_rh = alt.binding_select(options=[True, False, None], name="Real HETs")
selection_rh = alt.selection_single(fields=["realHET"], bind=selector_rh)

# Add Drop Down Selectors
vcTypes = ["INS", "DEL", "SNV", None]
illu_dd = alt.binding_select(options=vcTypes, name="Illumina VC Type")
illu_selection = alt.selection_single(fields=["illuminaType"], bind=illu_dd)
nano_dd = alt.binding_select(options=vcTypes, name="Nanopore VC Type")
nano_selection = alt.selection_single(fields=["nanoporeType"], bind=nano_dd)

# Combined Selection Interface
make_selector = (
    alt.Chart(make)
    .mark_rect()
    .encode(x="sample", color=color)
    .add_selection(selection_sm)
    .add_selection(selection_fp)
    .add_selection(selection_fn)
    .add_selection(selection_rh)
    .add_selection(selection_dh)
    .add_selection(illu_selection)
    .add_selection(nano_selection)
)

chart = (
    alt.Chart(df)
    .mark_rect()
    .encode(
        y="position:O",
        x="sample:N",
        color=alt.Color(
            "classification",
            scale=alt.Scale(
                domain=["Dropout", "Concordance", "Discordance"],
                range=["black", "blue", "orange"],
            ),
            legend=alt.Legend(orient="left"),
        ),
        tooltip=[
            "sample",
            "position",
            "reference",
            "illuminaType",
            "illuminaValue",
            "nanoporeType",
            "nanoporeValue",
            "nanoporeDropout",
            "illuminaDropout",
            "falseNegative",
            "falsePositive",
            "realHET",
            "detectedHET",
            "concordance",
            "illuminaPileup",
            "nanoporePileup",
        ],
    )
    .transform_filter(selection_sm)
    .transform_filter(selection_fn)
    .transform_filter(selection_fp)
    .transform_filter(selection_dh)
    .transform_filter(selection_rh)
    .transform_filter(illu_selection)
    .transform_filter(nano_selection)
    .interactive()
)

alt.vconcat(make_selector, chart, padding=256).save(snakemake.output[0])

with open(snakemake.output[0], "a") as outfile:
    outfile.write(
        " \
	<style>\
	form.vega-bindings {\
	  position: absolute;\
	  left: 0px;\
	  top: 0px;\
	}\
	</style>\
	"
    )
