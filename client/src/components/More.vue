<template>
    <div>
        <h1>ID: {{ $route.params.id }}</h1>
        <h4>Start: {{ currentCDS.start }}, Stop: {{ currentCDS.stop }}</h4>
        <br>
        <div align="center">
            <div class="alert alert-warning">
                <h3><strong>Status: </strong>{{ currentCDS.status }}</h3>
            </div>
        </div>
        <hr>
        <div id="docs" align="left">
        </div>
        <br>
        <div class="columns">
            <div class="column table">
                <div align="left" style="margin-left:37px; margin-right:50px;">
                    <h4 align="center">Coding Potential</h4>
                    <br>
                    <p>
                        The graph below displays the coding potential in each open reading frame (ORF). The first three plots represent direct sequences, and the latter three represent complementary sequences. See the table below the graph for more information.
                    </p>
                </div>
                <div id="direct-graph"></div>
                <div id="comp-graph"></div>
                <br><br>
                <div>
                    <table id="cp-table" class="table table-striped">
                        <caption>
                            This table displays the average coding potential per frame between the specified start position and {{ currentCDS.stop }} (the graph above is a visual representation of this data). Frames are counted 1-6 (direct 1-3 and complementary 4-6). Select a new start position, if appropriate. 
                        </caption>
                        <thead>
                            <div id="cp-head">
                                <tr>
                                    <th scope="col">Start Position</th>
                                    <th scope="col">Frame</th>
                                    <th scope="col">Average Coding Potential</th>
                                </tr>
                            </div>
                        </thead>
                        <tbody>
                            <div id="cp-body">
                                <tr>
                                    <td rowspan="0">
                                        {{ currentCDS.start }}
                                    </td>
                                    <td>1</td>
                                    <td>{{ probabilities.frame_1 }}</td>
                                </tr>
                                <tr>
                                    <td>2</td>
                                    <td>{{ probabilities.frame_2 }}</td>
                                </tr>
                                <tr>
                                    <td>3</td>
                                    <td>{{ probabilities.frame_3 }}</td>
                                </tr>
                                <tr>
                                    <td>4</td>
                                    <td>{{ probabilities.frame_4 }}</td>
                                </tr>
                                <tr>
                                    <td>5</td>
                                    <td>{{ probabilities.frame_5 }}</td>
                                </tr>
                                <tr>
                                    <td>6</td>
                                    <td>{{ probabilities.frame_6 }}</td>
                                </tr>
                            </div>
                        </tbody>
                    </table>
                    <br>
                </div>
            </div>
            <div class="column blast">
                <h4 align="center">BLAST</h4>
                <br>
                <h5>Instructions</h5>
                <p>
                    For this gene call, we will do a BLAST search on the sequence ranging from positions 
                    <strong>{{ currentCDS.start }}</strong> to <strong>{{ currentCDS.stop }}</strong>. 
                    You MUST stay on this page as you wait for results to appear.
                    Click the button below to begin the search.
                </p>
                <button class="btn btn-outline-primary" @click="runBLAST(currentCDS.id)">
                    BLAST
                    <b-spinner small v-if="showSpinner"><div class="divider"></div></b-spinner>
                </button>
                <br><br>
                <alert :message="message" v-if="showMessage"></alert>
                <br>
                <div id="blast-results" v-if="showBLAST">
                    <table id="blast-table" class="table table-hover">
                        <thead>
                            <tr>
                                <th id="hitDefHeader" scope="col">Description</th>
                                <th scope="col">E value</th>
                                <th scope="col">Query</th>
                                <th scope="col">Subject</th>
                            </tr>
                        </thead>
                        <tbody>
                            <tr v-for="(alignment, index) in blastResults" :key="index">
                                <td id="hitDef">{{ alignment.hit_def }}</td>
                                <td>{{ alignment.e_value }}</td>
                                <td>{{ alignment.query_start }} - {{ alignment.query_end }}</td>
                                <td>{{ alignment.sbjct_start }} - {{ alignment.sbjct_end }}</td>
                            </tr>
                        </tbody>
                    </table>
                </div>
                <div align="center">
                    <br><br>
                    <button class="btn btn-warning" @click="showModal">Select a function</button><br><br>
                    <span v-if="showFunction">Function: {{ newFunction }}</span><br>
                    <select id="functionSelect" v-model="newFunction" v-if="showFunction">
                        <option disabled value="">Please select one</option>
                        <option v-for="(alignment, index) in blastResults" :key="index">
                            {{ alignment.hit_def }}
                        </option>
                    </select>
                </div>
            </div>
        </div>
        <br>
        <button type="button" class="btn btn-success" :disabled="newFunction == ''" @click="editCDS">
            <strong>Update {{ $route.params.id }}</strong>
        </button>
        <br><br>
        <button type="button" class="btn btn-danger" @click="deleteCDS(currentCDS.id)">
            <strong>Delete {{ $route.params.id }}</strong>
        </button> 
        <br><br>
        <button class="btn btn-primary" @click="keepOriginal"><strong>Keep original data</strong></button>
    </div>
</template>

<script>
import axios from 'axios';
import Alert from './Alert.vue';
import {default as vegaEmbed} from 'vega-embed'

export default {
    data() {
        return {
            blastResults: [],
            showBLAST: false,
            message: '',
            showMessage: false,
            showSpinner: false,
            directGraph: [],
            compGraph: [],
            showFunction: false,
            showStart: false,
            newFunction: '',
            newStart: null,
            probabilities: {
                frame_1: 0,
                frame_2: 0,
                frame_3: 0,
                frame_4: 0,
                frame_5: 0,
                frame_6: 0,
            },
            currentCDS: {
                id: '',
                start: '',
                stop: '',
                strand: '',
                status: '',
            },
            updatedCDS: {
                id: '',
                start: '',
                stop: '',
                strand: '',
                function: '',
                status: '',
            },
        }
    },
    components: {
        alert: Alert,
    },
    methods: {
        getData(cdsID) {
            axios.get(`http://localhost:5000/annotate_data/more/${cdsID}`)
            .then(response => {
                this.currentCDS = response.data.cds;
                this.probabilities = response.data.probs;
                this.directGraph = response.data.direct_graph;
                this.compGraph = response.data.comp_graph;
                vegaEmbed('#direct-graph', this.directGraph, {actions: false})
                vegaEmbed('#comp-graph', this.compGraph, {actions: false})

            })
            .catch(error => {
                console.error(error);
            });
        },
        runBLAST(cdsID) {
            this.message = 'May take a few minutes. Please wait.'
            this.showMessage = true;
            this.showSpinner = true;
            axios.post(`http://localhost:5000/annotate_data/more/${cdsID}`)
            .then(response => {
                this.blastResults = response.data.blast;
                this.getData(this.$route.params.id);
                this.showBLAST = true;
                this.showMessage = true;
                this.message = 'BLAST complete!';
                this.showSpinner = false;
            })
            .catch(error => {
                console.error(error);
            });
        },
        showModal() {
            this.showFunction = true;
        },
        editCDS() {
            this.updatedCDS = this.currentCDS;
            this.updatedCDS.start = this.newStart;
            this.updatedCDS.function = this.newFunction;
            const payload = {
                id: this.updatedCDS.id,
                start: this.updatedCDS.start,
                stop: this.updatedCDS.stop,
                strand: this.updatedCDS.strand,
                function: this.updatedCDS.function,
                status: 'Pass',
            };
            this.updateCDS(payload, this.updatedCDS.id);
        },
        updateCDS(payload, cdsID) {
            axios.put(`http://localhost:5000/annotate_data/more/${cdsID}`, payload)
            .then(() => {
                this.$router.push("/annotate_data");
            })
            .catch(error => {
                console.error(error);
            });
        },
        deleteCDS(cdsID) {
            axios.delete(`http://localhost:5000/annotate_data/more/${cdsID}`)
            .then(() => {
                this.$router.push("/annotate_data");
            })
            .catch(error => {
                console.error(error);
                // this.getData();
            });
        },
        keepOriginal() {
            const payload = {
                id: this.currentCDS.id,
                start: this.currentCDS.start,
                stop: this.currentCDS.stop,
                strand: this.currentCDS.strand,
                function: 'None',
                status: 'Pass',
            };
            this.updateCDS(payload, this.currentCDS.id);
        },
    },
    created() {
        this.getData(this.$route.params.id);
    }
}
</script>

<style scoped>
.alert {
   width: 470px;
   height: 60px;
   text-align: center;
}

#statusMoreHeader {
    display: inline;
    border: 2px solid #ffbf00;
    border-radius: 15px;
    padding: 5px;
    color: #ffbf00;
    /* background-color: #d9eefc; */
}

#hitDefHeader {
    width: 50%;
}

#hitDef {
    width: 50%;
}

.divider {
    width:5px;
    height:auto;
    display:inline-block;
}

.columns {
    height: 100%;
    display: flex;
    /* vertical-align: stretch; */
}

#cp-table {
    width: 80%;
    margin-left: auto;
    margin-right: auto;
}

#cp-table tbody {
    overflow: auto;
    height:40.5vh;
    display: block;
}

#cp-head th {
    width:30vh;
}

#cp-body td {
    width:30vh;
}

caption {
  display: table-caption;
  caption-side: top;
}

.column.blast {
    flex-basis: 80%;
    width: 100vh;
    text-align: left;
}

#blast-results {
    overflow: auto;
    height: 211vh;
    flex-basis: 80%;
}
</style>