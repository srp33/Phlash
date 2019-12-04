<template>
    <div class="container">
        <h1>Database</h1>
        <h2>DNA Master</h2>
        <div id="dnamaster">
            <hr><br><br>
            <div align="left">
                <p>DOCS</p>
                <p>FIXME:</p>
                <ul>
                    <li>Delete button: are you sure you want to delete alert</li>
                    <li>When updating CDS, GET method doesn't resort start positions</li>
                </ul>
            </div>
            <br>
            <router-link to="/api/upload_genemark">
                <button class="btn btn-success"><strong>Next: Upload GeneMark file</strong></button>
            </router-link>
            <br><br>
            <alert :message="message" v-if="showMessage"></alert>
            <button class="btn btn-primary btn-sm" align="center" v-b-modal.cds-add-modal>Add a new CDS</button>
            <br><br>
            <table class="table table-hover" align="center">
                <thead>
                    <tr>
                        <th scope="col">ID</th>
                        <th scope="col">Start</th>
                        <th scope="col">Stop</th>
                        <th scope="col">Strand</th>
                        <th scope="col">Action</th>
                    </tr>
                </thead>
                <tbody>
                    <tr v-for="(curr, index) in dnamaster" :key="index">
                        <td>{{ curr.id }}</td>
                        <td>{{ curr.start }}</td>
                        <td>{{ curr.stop }}</td>
                        <td>{{ curr.strand }}</td>
                        <td>
                            <div class="btn-group" role="group">
                                <button type="button" class="btn btn-warning btn-sm" v-b-modal.cds-update-modal @click="editCDS(curr)">Update</button>
                                <div class="divider"></div>
                                <button type="button" class="btn btn-danger btn-sm" @click="removeCDS(curr.id)">Delete</button>
                            </div>
                        </td>
                    </tr>
                </tbody>
            </table>
        </div>
        <router-link to="/api/upload_genemark">
            <button class="btn btn-success"><strong>Next: Upload GeneMark file</strong></button>
        </router-link>
        <b-modal ref="addCDSModal" id="cds-add-modal" title="Add a new CDS" hide-footer>
            <b-form @submit="onSubmit" @reset="onReset" class="w-100" align="left">
                <b-form-group id = "form-id-group" label="ID:" label-for="form-id-input">
                    <b-form-input id="form-id-input" type="text" v-model="newCDS.id" required placeholder="Enter ID"></b-form-input>
                </b-form-group>
                <b-form-group id = "form-start-group" label="Start:" label-for="form-start-input">
                    <b-form-input id="form-start-input" type="text" v-model="newCDS.start" required placeholder="Enter start"></b-form-input>
                </b-form-group>
                <b-form-group id = "form-stop-group" label="Stop:" label-for="form-stop-input">
                    <b-form-input id="form-stop-input" type="text" v-model="newCDS.stop" required placeholder="Enter stop"></b-form-input>
                </b-form-group>
                <b-form-group id = "form-strand-group" label="Strand:" label-for="form-strand-input">
                    <b-form-input id="form-strand-input" type="text" v-model="newCDS.strand" required placeholder="Enter strand"></b-form-input>
                </b-form-group>
                <b-button type="submit" variant="primary">Submit</b-button>
                <b-button type="reset" variant="danger">Cancel</b-button>
            </b-form>
        </b-modal>
        <b-modal ref="updateCDSModal" id="cds-update-modal" title="Update CDS" hide-footer>
            <b-form @submit="onSubmitUpdate" @reset="onResetUpdate" class="w-100" align="left">
                <b-form-group id="form-start-update-group" label="Start:" label-for="form-start-update-input">
                    <b-form-input id="form-start-update-input" type="text" v-model="updatedCDS.start" required placeholder="Enter start"></b-form-input>
                </b-form-group>
                <b-form-group id="form-stop-update-group" label="Stop:" label-for="form-stop-update-input">
                    <b-form-input id="form-stop-update-input" type="text" v-model="updatedCDS.stop" required placeholder="Enter stop"></b-form-input>
                </b-form-group>
                <b-form-group id="form-strand-update-group" label="Strand:" label-for="form-strand-update-input">
                    <b-form-input id="form-strand-update-input" type="text" v-model="updatedCDS.strand" required placeholder="Enter strand"></b-form-input>
                </b-form-group>
                <b-button type="submit" variant="primary">Submit</b-button>
                <b-button type="reset" variant="danger">Cancel</b-button>
            </b-form>
        </b-modal>
    </div>
</template>

<script>
import axios from 'axios';
import Alert from './Alert.vue';

export default {
    data() {
        return {
            dnamaster: [],
            newCDS: {
                id: '',
                start: '',
                stop: '',
                strand: '',
                read: [],
            },
            updatedCDS: {
                id: '',
                start: '',
                stop: '',
                strand: '',
                read: [],
            },
            message: '',
            showMessage: false,
        };
    },
    components: {
        alert: Alert,
    },
    methods: {
        getData() {
            axios.get('http://localhost:5000/database')
            .then(response => {
                this.dnamaster = response.data.dnamaster;
            })
            .catch(error => {
                console.error(error);
            });
        },
        addData(payload) {
            axios.post("http://localhost:5000/database", payload)
            .then(() => {
                this.getData();
            })
            .catch(error => {
                console.log(error);
                this.getData();
            });
        },
        removeCDS(cdsID) {
            axios.delete(`http://localhost:5000/database/${cdsID}`)
            .then(() => {
                this.getData();
                this.message = 'CDS removed!';
                this.showMessage = true;
            })
            .catch(error => {
                console.error(error);
                this.getData();
            });
        },
        editCDS(cds) {
            this.updatedCDS.id = cds.id;
            this.updatedCDS.start = cds.start;
            this.updatedCDS.stop = cds.stop;
            this.updatedCDS.strand = cds.strand;
        },
        initForm() {
            this.newCDS.id = '';
            this.newCDS.start = '';
            this.newCDS.stop = '';
            this.newCDS.strand = '';
            this.updatedCDS.id = '';
            this.updatedCDS.start = '';
            this.updatedCDS.stop = '';
            this.updatedCDS.strand = '';
        },
        onSubmit(evt) {
            evt.preventDefault();
            this.$refs.addCDSModal.hide();
            let read = false;
            if (this.newCDS.read[0]) read = true;
            const payload = {
                id: this.newCDS.id,
                start: this.newCDS.start,
                stop: this.newCDS.stop,
                strand: this.newCDS.strand,
                read, // property shorthand
            };
            this.addData(payload);
            this.initForm();
            this.message = "CDS added!";
            this.showMessage = true;
        },
        onReset(evt) {
            evt.preventDefault();
            this.$refs.addCDSModal.hide();
            this.initForm();
            this.getData();
        },
        onSubmitUpdate(evt) {
            evt.preventDefault();
            this.$refs.updateCDSModal.hide();
            let read = false;
            if (this.updatedCDS.read[0]) read = true;
            const payload = {
                id: this.updatedCDS.id,
                start: this.updatedCDS.start,
                stop: this.updatedCDS.stop,
                strand: this.updatedCDS.strand,
                read, // property shorthand
            };
            this.updateCDS(payload, this.updatedCDS.id);
            this.message = 'CDS updated!';
            this.showMessage = true;
        },
        updateCDS(payload, cdsID) {
            axios.put(`http://localhost:5000/database/${cdsID}`, payload)
            .then(() => {
                this.getData();
            })
            .catch(error => {
                console.error(error);
                // this.getData();
            });
        },
        onResetUpdate(evt) {
            evt.preventDefault();
            this.$refs.updateCDSModal.hide();
            this.initForm();
            this.getData();
        },
    },
    created() {
        this.getData();
    },
};
</script>

<style scoped>
.divider {
    width:5px;
    height:auto;
    display:inline-block;
}
</style>