<template>
  <div class="navbar-wrapper">
    <nav
      class="navbar sticky-top navbar-expand-lg navbar-dark bg-dark"
      style="color: white; font-size: 1.5em"
    >
      <a class="navbar-brand" href="/">
        <img id="logo" src="/phlash/images/logonav.png" width="50" />
      </a>
      <div class="collapse navbar-collapse" id="navbarNav">
        <ul class="navbar-nav">
          <li class="nav-item dropdown">
            <a class="nav-link dropdown-toggle" href="#" id="navbarDropdown" role="button" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
              phlash
            </a>
            <div class="dropdown-menu" aria-labelledby="navbarDropdown">
              <a class="dropdown-item" href="#"><router-link to="/" style="color: black">home</router-link></a>
              <a v-if="settings" class="dropdown-item" href="#" @click="getSettings">settings</a>
              <a v-if="annotations && blast" class="dropdown-item" href="#" @click="showShare = true">share</a>
              <div v-if="logout" class="dropdown-divider"></div>
              <a v-if="logout" class="dropdown-item" style="margin-left:0px; text-align:left;" href="#"><GoogleLogin class="btn btn-dark btn-block" :params="params" :onSuccess="onSuccess" :logoutButton=true>logout</GoogleLogin></a>
            </div>
          </li>
          <li class="nav-item" v-if="upload">
            <router-link
              :to="{ name: 'Upload', params: { phageID: this.phageID } }"
              :event="upload ? 'click' : ''"
            >
              <a class="nav-link" href="#">upload</a>
            </router-link>
          </li>
          <li class="nav-item" v-if="blast">
            <router-link
              :to="{ name: 'Blast', params: { phageID: this.phageID } }"
              :event="blast ? 'click' : ''"
            >
              <a class="nav-link" href="#">blast</a>
            </router-link>
          </li>
          <li class="nav-item" v-if="annotations">
            <router-link
              :to="{ name: 'Annotations', params: { phageID: this.phageID } }"
              :event="annotations ? 'click' : ''"
            >
              <a class="nav-link" href="#">annotations</a>
            </router-link>
          </li>
          <li class="nav-item" v-if="annotations">
            <router-link
              :to="{ name: 'GeneMap', params: { phageID: this.phageID } }"
              :event="annotations ? 'click' : ''"
            >
              <a class="nav-link" href="#">genome map</a>
            </router-link>
          </li>
          <li class="nav-item" v-if="annotations">
            <router-link
              :to="{ name: 'GenBank', params: { phageID: this.phageID } }"
              :event="annotations ? 'click' : ''"
            >
              <a class="nav-link" href="#">genbank</a>
            </router-link>
          </li>
          <li class="nav-item">
            <router-link to="/contact"
              ><a class="nav-link" href="#">about</a></router-link
            >
          </li>
        </ul>
      </div>
      <span v-if="annotations && !blast" class="navbar-text">&#128065; VIEW ONLY &#128065;</span>
      <button
        class="navbar-toggler"
        type="button"
        data-toggle="collapse"
        data-target="#navbarNav"
        aria-controls="navbarNav"
        aria-expanded="false"
        aria-label="Toggle navigation"
      >
        <span class="navbar-toggler-icon"></span>
      </button>
    </nav>
    <b-modal
      class="text-size"
      v-model="showShare"
      ref="shareModal"
      id="share-modal"
      hide-footer
    >
      <template #modal-title>
        <div class="text-size">Share Annotations</div>
      </template>
      <b-form @submit="onShare" align="left">
        Enter the email of the user that you would like to share this phage's annotations with. 
        Please note that they will only be able to view your annotations and will not have permission to edit.
        <hr />
        <b-form-group
          label="Email:"
          label-size="lg"
          label-for="share-email"
        >
          <b-form-input
            class="form-input"
            id="share-email"
            type="text"
            v-model="shareEmail"
            required
            placeholder="jane_doe@gmail.com"
          ></b-form-input>
        </b-form-group>
        <hr />
        <b-button type="submit" class="mt-3" block style="margin-top: 0em">
          <strong>Submit</strong>
        </b-button>
      </b-form>
    </b-modal>
    <b-modal
      class="text-size"
      v-model="showSettings"
      ref="settingsModal"
      id="settings-modal"
      hide-footer
    >
      <template #modal-title>
        <div class="text-size">Settings</div>
      </template>
      <b-form @submit="onSubmit" align="left">
        <button class="btn btn-secondary btn-block" style="margin-top:1em; margin-bottom:.5em;" type="button" data-toggle="collapse" data-target="#AnnotationSoftware" aria-expanded="false" aria-controls="AnnotationSoftware">
          Annotation Software Settings
        </button>
        <div class="collapse" id="AnnotationSoftware">
          Click 'Submit' and remove and re-upload the fasta file for the new settings to take effect.
          <hr />
          <b-form-checkbox
            id="glimmer-input"
            type="checkbox"
            v-model="glimmer"
          > Glimmer3
          </b-form-checkbox>
          <b-form-checkbox
            id="genemark-input"
            type="checkbox"
            v-model="genemark"
          > GeneMarkS
          </b-form-checkbox>
          <b-form-checkbox
            id="prodigal-input"
            type="checkbox"
            v-model="prodigal"
          > Prodigal
          </b-form-checkbox>
          <b-form-checkbox
            id="phanotate-input"
            type="checkbox"
            v-model="phanotate"
          > PHANOTATE
          </b-form-checkbox>
          <b-form-checkbox
            id="aragorn-input"
            type="checkbox"
            v-model="aragorn"
          > ARAGORN
          </b-form-checkbox>
        </div>
        <button class="btn btn-secondary btn-block" style="margin-top:1em; margin-bottom:.5em;" type="button" data-toggle="collapse" data-target="#SearchRange" aria-expanded="false" aria-controls="SearchRange">
          Search Range Settings
        </button>
        <div class="collapse" id="SearchRange">
          Click 'Submit' and refresh the page for the new settings to take effect.
          <hr />
          <b-form-group
            label="Start Search Range (back):"
            label-size="lg"
            label-for="search-back-input"
          >
            This number (defaulted at 300) represents how many base pairs back
            from the current start position alternate start codons will be
            searched for.
            <b-form-input
              class="form-input"
              id="search-back-input"
              type="number"
              v-model="backStartRange"
              required
              placeholder="Enter start search range"
            ></b-form-input>
          </b-form-group>
          <b-form-group
            label="Start Search Range (forward):"
            label-size="lg"
            label-for="search-forward-input"
          >
            This number (defaulted at 100) represents how many base pairs forward
            from the current start position alternate start codons will be
            searched for.
            <b-form-input
              class="form-input"
              id="search-forward-input"
              type="number"
              v-model="forwardStartRange"
              required
              placeholder="Enter start search range"
            ></b-form-input>
          </b-form-group>
        </div>
        <button class="btn btn-secondary btn-block" style="margin-top:1em; margin-bottom:.5em;" type="button" data-toggle="collapse" data-target="#Status" aria-expanded="false" aria-controls="Status">
          Status Settings
        </button>
        <div class="collapse" id="Status">
          Click 'Submit' and refresh the page for the new settings to take effect.
          <hr />
          <b-form-group label="Gap:" label-size="lg" label-for="gap-input">
            This number (defaulted at 10) represents the maximum acceptable number
            of base pairs between two adjacent genes. Gaps greater than this
            number will be flagged.
            <b-form-input
              class="form-input"
              id="gap-input"
              type="number"
              v-model="gap"
              required
              placeholder="Enter acceptable gap length"
            ></b-form-input>
          </b-form-group>
          <b-form-group
            label="Overlap:"
            label-size="lg"
            label-for="overlap-input"
          >
            This number (defaulted at 10) represents the maximum acceptable number
            of base pairs two adjacent genes overlap. Overlaps greater than this
            number will be flagged.
            <b-form-input
              class="form-input"
              id="overlap-input"
              type="number"
              v-model="overlap"
              required
              placeholder="Enter acceptable overlap length."
            ></b-form-input>
          </b-form-group>
          <b-form-group
            label="Opposite Strand Gap:"
            label-size="lg"
            label-for="opposite-gap-input"
          >
            This number (defaulted at 50) represents the minimum acceptable number
            of base pairs between two adjacent genes on different strands. Gaps
            shorter than this number will be flagged.
            <b-form-input
              class="form-input"
              id="opposite-gap-input"
              type="number"
              v-model="oppositeGap"
              required
              placeholder="Enter acceptable gap length"
            ></b-form-input>
          </b-form-group>
          <b-form-group
            label="Minimum Gene Length:"
            label-size="lg"
            label-for="short-input"
          >
            This number (defaulted at 200) represents the minimum acceptable
            number of base pairs in a gene. Genes shorter than this number will be
            flagged.
            <b-form-input
              class="form-input"
              id="short-input"
              type="number"
              v-model="short"
              required
              placeholder="Enter acceptable minimum gene length"
            ></b-form-input>
          </b-form-group>
        </div>
        <hr />
        <b-button type="submit" class="mt-3" block style="margin-top: 0em">
          <strong>Submit</strong>
        </b-button>
      </b-form>
    </b-modal>
    <b-toast id="share-status" variant="primary">
      <template #toast-title>
        <strong class="text-size"> {{statusTitle}} </strong>
      </template>
      <div class="text-size">{{ statusMessage }}</div>
    </b-toast>
  </div>
</template>

<script>
import axios from 'axios';
import GoogleLogin from 'vue-google-login';

export default {
  name: 'Navbar',
  components: {
    GoogleLogin,
  },
  props: {
    upload: Boolean,
    blast: Boolean,
    annotations: Boolean,
    geneMap: Boolean,
    settings: Boolean,
    phageID: String,
    logout: Boolean,
  },

  data() {
    return {
      activeCLass: 'active',
      showSettings: false,
      prodigal: null,
      glimmer: null,
      genemark: null,
      aragorn: null,
      phanotate: null,
      gap: null,
      overlap: null,
      oppositeGap: null,
      backStartRange: null,
      forwardStartRange: null,
      short: null,
      showShare: false,
      shareEmail: null,
      statusMessage: "",
      statusTitle: "",
    };
  },

  methods: {
    getSettings() {
      (this.showSettings = true),
        axios
          .options(
            process.env.VUE_APP_BASE_URL +
              `/annotations/${this.$route.params.phageID}/none`
          )
          .then((response) => {
            this.prodigal = response.data.prodigal;
            this.glimmer = response.data.glimmer;
            this.genemark = response.data.genemark;
            this.aragorn = response.data.aragorn;
            this.phanotate = response.data.phanotate;
            this.backStartRange = response.data.back_left_range;
            this.forwardStartRange = response.data.forward_left_range;
            this.gap = response.data.gap;
            this.overlap = response.data.overlap;
            this.oppositeGap = response.data.opposite_gap;
            this.short = response.data.short;
          });
    },

    onSubmit(evt) {
      evt.preventDefault();
      this.$refs.settingsModal.hide();
      var payload =
        this.gap +
        ',' +
        this.overlap +
        ',' +
        this.oppositeGap +
        ',' +
        this.backStartRange +
        ',' +
        this.forwardStartRange +
        ',' +
        this.short +
        ',' +
        this.prodigal.toString() +
        ',' + 
        this.glimmer.toString() +
        ',' + 
        this.genemark.toString() +
        ',' + 
        this.aragorn.toString() +
        ',' + 
        this.phanotate.toString();
      axios
        .options(
          process.env.VUE_APP_BASE_URL +
            `/annotations/${this.$route.params.phageID}/${payload}`
        )
        .then((response) => {
          console.log(response.data);
        })
        .catch((error) => {
          console.error(error);
        });
    },

    onShare(evt) {
      evt.preventDefault();
      this.$refs.shareModal.hide();
      axios
        .post(
          process.env.VUE_APP_BASE_URL +
            `/annotations/${this.$route.params.phageID}/${this.shareEmail}`
        )
        .then((response) => {
          this.statusTitle = "STATUS";
          this.statusMessage = response.data;
          this.$bvToast.show('share-status');
        })
        .catch((error) => {
          console.error(error);
        });
    },

    onSuccess() {
      window.location.reload();
    },
  },
};
</script>

<style scoped>
.active {
  opacity: 1;
  visibility: visible;
  border-left-color: #4db6ac;
  margin: 0.25em;
  transition: all 0.25s;
}

.bg-light {
  background-color: #e3f2fd;
}

.text-size {
  font-size: 1.2em;
}

.form-input {
  height: 2em; 
  font-size: 15pt;
}

</style>