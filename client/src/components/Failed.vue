<template>
    <div>
        <h1>ID: {{ $route.params.id }}</h1>
        <h4>Start: {{ currentCDS.start }}, Stop: {{ currentCDS.stop }}</h4>
        <br>
        <div align="center">
            <div class="alert alert-danger" id="status-header">
                <h3><strong>Status: </strong>{{ currentCDS.status }}</h3>
            </div>
        </div>
        <hr>
        <div id="docs">
            <h5>Other possible start position(s):</h5>
            <h5 v-for="(curr, index) in startOptions" :key=index>
                {{ curr.start_position }}
            </h5>
        </div>
        <br>
        <div class="columns">
            <div class="column cp">
                <div align="left" style="margin-left:37px; margin-right:50px;">
                    <h4 align="center">Coding Potential</h4>
                    <br>
                    <p>
                        The graph below displays the coding potential in each open reading frame (ORF). The first three plots represent direct sequences, and the latter three represent complementary sequences. See the table below the graph for more information.
                    </p>
                    <h5>Instructions</h5>
                    <p style="color:#0066CC">
                        Hover your mouse over the graph to see the exact base positions. Then, scroll down to the table below and <strong>select a new start position, if appropriate</strong>.
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
                            <div id="cp-body" v-for="(curr, index) in startOptions" :key="index">
                                <tr>
                                    <td rowspan="0">
                                        {{ curr.start_position }}<br />
                                        <button class="btn btn-warning btn-sm" @click="setStart(curr.start_position)">Select</button>
                                    </td>
                                    <td>1</td>
                                    <td>{{ curr.frame_1 }}</td>
                                </tr>
                                <tr>
                                    <td>2</td>
                                    <td>{{ curr.frame_2 }}</td>
                                </tr>
                                <tr>
                                    <td>3</td>
                                    <td>{{ curr.frame_3 }}</td>
                                </tr>
                                <tr>
                                    <td>4</td>
                                    <td>{{ curr.frame_4 }}</td>
                                </tr>
                                <tr>
                                    <td>5</td>
                                    <td>{{ curr.frame_5 }}</td>
                                </tr>
                                <tr>
                                    <td>6</td>
                                    <td>{{ curr.frame_6 }}</td>
                                </tr>
                            </div>
                        </tbody>
                    </table>
                    <br>
                    <h5><strong>Start position selected: </strong>{{ newStart }}</h5>
                    <br>
                </div>
            </div>
            <div class="column blast">
                <div align="left">
                    <h4 align="center">BLAST</h4>
                    <br>
                    <p>
                        For this gene call, we will do a BLAST search on the sequence ranging from positions 
                        <strong>{{ lowestStart }}</strong> to <strong>{{ currentCDS.stop }}</strong>.
                    </p>
                    <h5>Instructions</h5>
                    <p style="color:#0066CC">
                        Click the button below to begin the search. You <strong>must</strong> stay on this page as you wait for results to appear. Once results appear, <strong>select the most appropriate function for this gene call</strong>.
                    </p>
                    <button class="btn btn-outline-primary" @click="runBLAST(currentCDS.id)" :disabled="disableBlAST">
                        BLAST
                        <b-spinner small v-if="showSpinner"><div class="divider"></div></b-spinner>
                    </button>
                </div>
                <br><br>
                <alert :message="message" v-if="showMessage"></alert>
                <br>
                <div v-if="blastResults.length">
                    <div id="blast-results" v-if="showBLAST">
                        <table id="blast-table" class="table table-hover">
                            <thead align="center">
                                <tr>
                                    <th scope="col" style="width:33%;">Description</th>
                                    <th scope="col">E-value</th>
                                    <th scope="col">Query</th>
                                    <th scope="col">Subject</th>
                                    <th scope="col">Action</th>
                                </tr>
                            </thead>
                            <tbody>
                                <tr v-for="(alignment, index) in blastResults" :key="index">
                                    <td align="left" style="width:33%;">{{ alignment.hit_def }}</td>
                                    <td align="center" style="width:15%;">{{ alignment.e_value }}</td>
                                    <td align="center">{{ alignment.query_start }} - {{ alignment.query_end }}</td>
                                    <td align="center">{{ alignment.sbjct_start }} - {{ alignment.sbjct_end }}</td>
                                    <td align="center">
                                        <button class="btn btn-warning btn-sm" @click="setFunction(alignment.hit_def)">Select</button>
                                    </td>
                                </tr>
                            </tbody>
                        </table>
                    </div>
                    <br>
                    <h5 align="center" v-if="showBLAST"><strong>Function selected: </strong>{{ newFunction }}</h5>
                </div>
                <div v-else>
                    <div class="alert alert-warning" v-if="showBLAST">
                        No results with E-value lower than threshold, <strong>1e-35</strong>.
                    </div>
                </div>
            </div>
        </div>
        <br>
        <button type="button" class="btn btn-success" :disabled="newFunction == '' || newStart == null" @click="editCDS">
            <strong>Update {{ $route.params.id }}</strong>
        </button>
        <br><br>
        <button class="btn btn-primary" @click="keepOriginal"><strong>Keep original data</strong></button>
    </div>
</template>

<script>
import axios from 'axios';
import Alert from './Alert.vue';
import {default as vegaEmbed} from 'vega-embed'
import VueApexCharts from 'vue-apexcharts';

export default {
    data () {
        return {
            startOptions: [],
            lowestStart: 0,
            currentCDS: {
                id: '',
                start: '',
                stop: '',
                strand: '',
                function: '',
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
            disableBlAST: true,
        }
    },
    components: {
        alert: Alert,
    },
    methods: {
        getData(cdsID) {
            axios.get(`http://localhost:5000/annotate/failed/${cdsID}`)
            .then(response => {
                this.currentCDS = response.data.cds;
                this.startOptions = response.data.start_options;
                this.lowestStart = response.data.lowest_start;
                this.directGraph = response.data.direct_graph;
                this.compGraph = response.data.comp_graph;
                vegaEmbed('#direct-graph', this.directGraph, {actions: false})
                vegaEmbed('#comp-graph', this.compGraph, {actions: false})
                this.disableBlAST = false;
            })
            .catch(error => {
                console.error(error);
            });
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
            axios.put(`http://localhost:5000/annotate/failed/${cdsID}`, payload)
            .then(() => {
                this.$router.push("/annotate_data");
            })
            .catch(error => {
                console.error(error);
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
        runBLAST(cdsID) {
            this.message = 'May take a few minutes. Please wait.'
            this.showMessage = true;
            this.showSpinner = true;
            axios.post(`http://localhost:5000/annotate/failed/${cdsID}`)
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
        setFunction(funct) {
            this.newFunction = funct;
        },
        setStart(strt) {
            this.newStart = strt;
        },
    },
    created() {
        this.getData(this.$route.params.id);
    }
}
</script>

<style scoped>
#status-header {
   width: 190px;
   height: 60px;
   text-align: center;
}

#docs {
    margin-left: 150px;
    margin-right: 150px;
}

.columns {
    height: 100%;
    display: flex;
}

.column.cp {
    margin-right: 2px;
    text-align: center;
}

.column.blast {
    flex-basis: 110%;
    width: 100vh; 
    text-align: left;
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
    width:26.3vh;
}

#cp-body td{
    width:30vh;
}

caption {
  display: table-caption;
  caption-side: top;
  /* text-align: center; */
}

#blast-table thead tr {
    display: block;
} 

#blast-table tbody {
    overflow: auto;
    height: 203vh;
    display: block;
}

#blast-table th {
    width: 98px;
}

#blast-table td {
    width: 96px;
}

.divider {
    width:5px;
    height:auto;
    display:inline-block;
}
</style>