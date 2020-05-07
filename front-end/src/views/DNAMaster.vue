<template>
  <div class="wrapper">
    <Navbar :upload="navUpload" :dnamaster="navDNAMaster" :blast="navBlast" :annotations="navAnnotations" />
    <div class="container">
      <h1>DNA Master</h1>
      <div class="alert alert-primary">
        <p><strong>Instructions</strong></p>
        <p>
          Below is the data from your DNA Master GenBank file. If you would like
          to modify any information <em>before beginning the annotation
          process</em>, please <b>add</b>, <b>update</b>, or <b>delete</b> the
          appropriate coding sequences.
        </p>
        <p>Click `Next` when you are ready to continue.</p>
        <div class="nav-btns-wrapper">
          <router-link :to="{ name: 'Upload', params: {phageID: $route.params.phageID} }">
            <button class="btn btn-light btn-nav">
              <svg class="bi bi-arrow-left" width="1em" height="1em" viewBox="0 0 16 16" fill="currentColor" xmlns="http://www.w3.org/2000/svg">
                <path fill-rule="evenodd" d="M5.854 4.646a.5.5 0 010 .708L3.207 8l2.647 2.646a.5.5 0 01-.708.708l-3-3a.5.5 0 010-.708l3-3a.5.5 0 01.708 0z" clip-rule="evenodd"/>
                <path fill-rule="evenodd" d="M2.5 8a.5.5 0 01.5-.5h10.5a.5.5 0 010 1H3a.5.5 0 01-.5-.5z" clip-rule="evenodd"/>
              </svg>
              <strong>Back</strong>
            </button>
          </router-link>
          <router-link :to="{ name: 'Blast', params: {phageID: $route.params.phageID} }">
            <button class="btn btn-light btn-nav">
              <strong>Next</strong>
              <svg class="bi bi-arrow-right" width="1em" height="1em" viewBox="0 0 16 16" fill="currentColor" xmlns="http://www.w3.org/2000/svg">
                <path fill-rule="evenodd" d="M10.146 4.646a.5.5 0 01.708 0l3 3a.5.5 0 010 .708l-3 3a.5.5 0 01-.708-.708L12.793 8l-2.647-2.646a.5.5 0 010-.708z" clip-rule="evenodd"/>
                <path fill-rule="evenodd" d="M2 8a.5.5 0 01.5-.5H13a.5.5 0 010 1H2.5A.5.5 0 012 8z" clip-rule="evenodd"/>
              </svg>
            </button>
          </router-link>
        </div>
      </div>

      <div class="alert alert-success" id="success-alert" role="alert" v-if="showMessage">
        {{ message }}
      </div>

      <div id="dnamaster">
        <div class="table-responsive">
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
                    <button class="btn btn-dark btn-sm" id="update-btn" v-b-modal.update-modal @click="editCDS(curr)">
                      <strong>Update</strong>
                    </button>
                    <button class="btn btn-dark btn-sm" id="delete-btn" v-b-modal.delete-modal @click="deleteCDS(curr)">
                      <strong>Delete</strong>
                    </button>
                  </div>
                </td>
              </tr>
            </tbody>
          </table>
        </div>
      </div>
      
      <div class="nav-btns-wrapper">
        <router-link :to="{ name: 'Upload', params: {phageID: $route.params.phageID} }">
          <button class="btn btn-light btn-nav">
            <svg class="bi bi-arrow-left" width="1em" height="1em" viewBox="0 0 16 16" fill="currentColor" xmlns="http://www.w3.org/2000/svg">
              <path fill-rule="evenodd" d="M5.854 4.646a.5.5 0 010 .708L3.207 8l2.647 2.646a.5.5 0 01-.708.708l-3-3a.5.5 0 010-.708l3-3a.5.5 0 01.708 0z" clip-rule="evenodd"/>
              <path fill-rule="evenodd" d="M2.5 8a.5.5 0 01.5-.5h10.5a.5.5 0 010 1H3a.5.5 0 01-.5-.5z" clip-rule="evenodd"/>
            </svg>
            <strong>Back</strong>
          </button>
        </router-link>
        <router-link :to="{ name: 'Blast', params: {phageID: $route.params.phageID} }">
          <button class="btn btn-light btn-nav">
            <strong>Next</strong>
            <svg class="bi bi-arrow-right" width="1em" height="1em" viewBox="0 0 16 16" fill="currentColor" xmlns="http://www.w3.org/2000/svg">
              <path fill-rule="evenodd" d="M10.146 4.646a.5.5 0 01.708 0l3 3a.5.5 0 010 .708l-3 3a.5.5 0 01-.708-.708L12.793 8l-2.647-2.646a.5.5 0 010-.708z" clip-rule="evenodd"/>
              <path fill-rule="evenodd" d="M2 8a.5.5 0 01.5-.5H13a.5.5 0 010 1H2.5A.5.5 0 012 8z" clip-rule="evenodd"/>
            </svg>
          </button>
        </router-link>
      </div>

      <b-modal ref="updateCDSModal" id="update-modal" title="Update CDS" hide-footer>
        <b-form @submit="onSubmitUpdate" @reset="onResetUpdate" align="left">
          <b-form-group label="Start:" label-for="update-start-input">
            <b-form-input
              id="update-start-input"
              type="number"
              v-model="updatedCDS.start"
              required
              placeholder="Enter start position"
            ></b-form-input>
          </b-form-group>
          <b-form-group label="Stop:" label-for="update-stop-input">
            <b-form-input
              id="update-stop-input"
              type="number"
              v-model="updatedCDS.stop"
              required
              placeholder="Enter stop position"
            ></b-form-input>
          </b-form-group>
          <b-form-group label="Strand:">
            <b-form-select v-model="updatedCDS.strand" :options="strandOptions"></b-form-select>
          </b-form-group>
          <hr />
          <b-button class="modal-btn" type="submit" variant="success">
            <strong>Submit</strong>
          </b-button>
          <b-button class="modal-btn" type="reset" variant="danger">
            <strong>Cancel</strong>
          </b-button>
        </b-form>
      </b-modal>
      <b-modal ref="deleteCDSModal" id="delete-modal" title="Warning" hide-footer>
        <b-form @submit="onSubmitDelete" @reset="onResetDelete" align="left">
          <p>Are you sure you want to delete this CDS?</p>
          <ul>
            <li><strong>ID:</strong> {{ deletedCDS.id }}</li>
            <li><strong>Start:</strong> {{ deletedCDS.start }}</li>
            <li><strong>Stop:</strong> {{ deletedCDS.stop }}</li>
            <li><strong>Strand:</strong> {{ deletedCDS.strand }}</li>
          </ul>
          <hr />
          <b-button class="modal-btn" type="submit" variant="success">
            <strong>Yes, delete it</strong>
          </b-button>
          <b-button class="modal-btn" type="reset" variant="danger">
            <strong>No, cancel</strong>
          </b-button>
        </b-form>
      </b-modal>
    </div>
  </div>
</template>

<script>
import axios from "axios";
import Navbar from "../components/Navbar.vue";

export default {
  name: "DNAMaster",
  components: {
    Navbar
  },
  data() {
    return {
      dnamaster: [],
      updatedCDS: {
        id: "",
        start: "",
        stop: "",
        strand: null,
        read: []
      },
      deletedCDS: {
        id: "",
        start: "",
        stop: "",
        strand: ""
      },
      strandOptions: [
        { value: null, text: "Please select a direction" },
        { value: "+", text: "+ (Direct)" },
        { value: "-", text: "- (Complementary)" }
      ],
      message: "",
      showMessage: false
    };
  },
  created() {
    this.getData();
  },
  computed: {
    navUpload: function() {
      return true;
    },
    navDNAMaster: function() {
      return true;
    },
    navBlast: function() {
      return true;
    },
    navAnnotations: function() {
      return false;
    },
  },
  methods: {
    getData() {
      axios.get(`http://localhost:5000/api/dnamaster/${this.$route.params.phageID}`)
        .then(response => {
          this.dnamaster = response.data.dnamaster;
        })
        .catch(error => {
          console.error(error);
        });
    },
    deleteData(cdsID) {
      axios.delete(`http://localhost:5000/api/dnamaster/${this.$route.params.phageID}/${cdsID}`)
        .then(response => {
          this.message = response.data.message;
          this.showMessage = true;
          this.getData();
        })
        .catch(error => {
          console.error(error);
        });
    },
    updateData(payload, cdsID) {
      axios.put(`http://localhost:5000/api/dnamaster/${this.$route.params.phageID}/${cdsID}`,
          payload
        )
        .then(response => {
          this.message = response.data.message;
          this.showMessage = true;
          this.getData();
        })
        .catch(error => {
          console.error(error);
        });
    },
    editCDS(cds) {
      this.updatedCDS.id = cds.id;
      this.updatedCDS.start = cds.start;
      this.updatedCDS.stop = cds.stop;
      this.updatedCDS.strand = cds.strand;
    },
    deleteCDS(cds) {
      this.deletedCDS.id = cds.id;
      this.deletedCDS.start = cds.start;
      this.deletedCDS.stop = cds.stop;
      this.deletedCDS.strand = cds.strand;
    },
    initForm() {
      this.updatedCDS.id = "";
      this.updatedCDS.start = "";
      this.updatedCDS.stop = "";
      this.updatedCDS.strand = "";
      this.deletedCDS.id = "";
      this.deletedCDS.start = "";
      this.deletedCDS.stop = "";
      this.deletedCDS.strand = "";
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
        read // property shorthand
      };
      this.updateData(payload, this.updatedCDS.id);
    },
    onResetUpdate(evt) {
      evt.preventDefault();
      this.$refs.updateCDSModal.hide();
      this.initForm();
      this.getData();
    },
    onSubmitDelete(evt) {
      evt.preventDefault();
      this.$refs.deleteCDSModal.hide();
      let read = false;
      this.deleteData(this.deletedCDS.id);
    },
    onResetDelete(evt) {
      evt.preventDefault();
      this.$refs.deleteCDSModal.hide();
      this.initForm();
      this.getData();
    }
  }
};
</script>

<style scoped>
.wrapper {
  margin: 0;
}

h1 {
  margin: 40px auto;
}

.alert-primary {
  text-align: left;
  margin: 40px auto;
}

.nav-btns-wrapper {
  text-align: center;
}

.btn-nav {
  margin: 5px;
}

.bi-arrow-left {
  margin-right: 5px;
  margin-left: 0px;
}

.bi-arrow-right {
  margin-right: 0px;
  margin-left: 5px;
}

#update-btn,
#delete-btn {
  margin: auto 3px;
}

/* ----- Modals ----- */
.modal-btn {
  margin: 3px;
}

.strand-dropdown {
  width: 100%;
  text-align: left;
  border: none;
  background: none;
}
</style>