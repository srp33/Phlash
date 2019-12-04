<template>
    <div id="blastForFunction">
        <h1>ID: {{ $route.params.id }}</h1>
        <h4>Start: {{ currentCDS.start }}, Stop: {{ currentCDS.stop }}</h4>
        <br>
        <div align="center">
            <div class="alert alert-success">
                <h3><strong>Status: </strong>{{ currentCDS.status }}</h3>
            </div>
        </div>
        <hr>
        <div align="left">
            <h3 align="center">BLAST for Function</h3><br>
            <p>
                For this gene call, we will do a BLAST search on the sequence ranging from positions 
                <strong>{{ lowestStart }}</strong> to <strong>{{ currentCDS.stop }}</strong>.
            </p>
            <h5>Instructions</h5>
            <p style="color:#0066CC">
                Click the button below to begin the search. You <strong>must</strong> stay on this page as you wait for results to appear. Once results appear, <strong>select the most appropriate function for this gene call</strong>.
            </p>
            <button class="btn btn-outline-primary" @click="runBLAST(currentCDS.id)">
                BLAST 
                <b-spinner small v-if="showSpinner">
                    <div class="divider"></div>
                </b-spinner>
            </button>
            <br><br>
            <alert :message="message" v-if="showMessage"></alert>
            <br>
        </div>
        <div id="blast-results" v-if="showBLAST">
            <br><br>
            <table id="blast-table" class="table table-hover">
                <thead>
                    <tr>
                        <th id="hitDefHeader" scope="col">Function</th>
                        <th scope="col">E value</th>
                        <th scope="col">Query</th>
                        <th scope="col">Subject</th>
                        <th scope="col">Action</th>
                    </tr>
                </thead>
                <tbody>
                    <tr v-for="(alignment, index) in blastResults" :key="index">
                        <td id="hitDef">{{ alignment.hit_def }}</td>
                        <td>{{ alignment.e_value }}</td>
                        <td>{{ alignment.query_start }} - {{ alignment.query_end }}</td>
                        <td>{{ alignment.sbjct_start }} - {{ alignment.sbjct_end }}</td>
                        <td>
                            <button class="btn btn-warning btn-sm" @click="setFunction(alignment.hit_def)">Select</button>
                        </td>
                    </tr>
                </tbody>
            </table>
            <br><br>
            <h4 v-if="showBLAST"><strong>Selection: </strong>{{ newFunction }}</h4>
        </div>
        <!-- <div align="center">
            <br><br>
            <button class="btn btn-warning" @click="showModal" v-if="showBLAST">Select a function</button><br><br>
            <span v-if="showFunction"><strong>Function: </strong>"{{ newFunction }}"</span><br>
            <select id="functionSelect" v-model="newFunction" v-if="showFunction">
                <option disabled value="">Please select one</option>
                <option v-for="(alignment, index) in blastResults" :key="index">
                    {{ alignment.hit_def }}
                </option>
            </select>
        </div> -->
        <br>
        <button type="button" class="btn btn-success" @click="editCDS" :disabled="newFunction == ''">
            <strong>Update {{ $route.params.id }}</strong>
        </button>
        <br><br>
        <button class="btn btn-primary" @click="keepOriginal"><strong>Go back</strong></button>
    </div>
</template>

<script>
import axios from 'axios';
import Alert from './Alert.vue';

export default {
    data () {
        return {
            blastResults: [],
            showBLAST: false,
            message: '',
            showMessage: false,
            showSpinner: false,
            showFunction: false,
            newFunction: '',
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
        }
    },
    components: {
        alert: Alert,
    },
    methods: {
        getData(cdsID) {
            axios.get(`http://localhost:5000/blast/${cdsID}`)
            .then(response => {
                this.currentCDS = response.data.cds;
            })
            .catch(error => {
                console.error(error);
            });
        },
        runBLAST(cdsID) {
            this.message = 'May take a few minutes. Please wait.'
            this.showMessage = true;
            this.showSpinner = true;
            axios.post(`http://localhost:5000/blast/${cdsID}`)
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
        editCDS() {
            this.updatedCDS = this.currentCDS;
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
            axios.put(`http://localhost:5000/blast/${cdsID}`, payload)
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
                function: this.currentCDS.function,
                status: this.currentCDS.status,
            };
            this.updateCDS(payload, this.currentCDS.id);
        },
    },
    created() {
        this.getData(this.$route.params.id);
    },
}
</script>

<style scoped>
#blastForFunction {
    margin-left: 100px;
    margin-right: 100px;
}

.alert {
   width: 200px;
   height: 60px;
   text-align: center;
}

#hitDefHeader {
    width: 40%;
}

#hitDef {
    width: 40%;
    text-align: left;
}

#blast-table thead tr {
    display: block;
} 

#blast-table tbody {
    overflow: auto;
    height: 70vh;
    display: block;
}

th, td {
    width: 166px;
}

.divider {
    width:5px;
    height:auto;
    display:inline-block;
}
</style>