<template>
  <div class="navbar-wrapper">
    <nav class="navbar sticky-top navbar-expand-lg navbar-light bg-light">
      <a class="navbar-brand" href="/">Phlash</a>
      <button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbarNav" aria-controls="navbarNav" aria-expanded="false" aria-label="Toggle navigation">
        <span class="navbar-toggler-icon"></span>
      </button>
      <div class="collapse navbar-collapse" id="navbarNav">
        <ul class="navbar-nav">
          <li class="nav-item">
            <router-link to="/"><a class="nav-link" href="#">home</a></router-link>
          </li>
          <li class="nav-item">
            <router-link :to="{ name: 'Upload', params: {phageID: this.phageID} }"
                         :event="upload ? 'click' : ''">
              <a class="nav-link" href="#">upload</a>
            </router-link>
          </li>
          <!-- <li class="nav-item">
            <router-link :to="{ name: 'DNAMaster', params: {phageID: this.phageID} }"
                         :event="dnamaster ? 'click' : ''">
              <a class="nav-link" href="#">dna master</a>
            </router-link>
          </li> -->
          <li class="nav-item">
            <router-link :to="{ name: 'Blast', params: {phageID: this.phageID} }"
                         :event="geneMap ? 'click' : ''">
              <a class="nav-link" href="#">blast</a>
            </router-link>
          </li>
          <li class="nav-item">
            <router-link :to="{ name: 'Annotations', params: {phageID: this.phageID} }"
                         :event="annotations ? 'click' : ''">
              <a class="nav-link" href="#">annotations</a>
            </router-link>
          </li>
          <li class="nav-item">
            <router-link :to="{ name: 'GeneMap', params: {phageID: this.phageID} }"
                         :event="blast ? 'click' : ''">
              <a class="nav-link" href="#">gene map</a>
            </router-link>
          </li>
          <li class="nav-item" v-if="settings">
            <a class="nav-link" href="#" @click="getSettings">settings</a>
          </li>
        </ul>
      </div>
    </nav>
      <b-modal v-model="showSettings" ref="settingsModal" id="settings-modal" title="Settings" hide-footer>
        Refresh the page for the new settings to take effect.
        <hr />
        <b-form @submit="onSubmit" align="left">
          <b-form-group label="Start Search Range (back):" label-size="lg" label-for="search-back-input">
            This number (defaulted at 300) represents how many base pairs back from the current start position 
            alternate start codons will be searched for.
            <b-form-input
              id="search-back-input"
              type="number"
              v-model="backStartRange"
              required
              placeholder="Enter start search range"
            ></b-form-input>
          </b-form-group>
          <b-form-group label="Start Search Range (forward):" label-size="lg" label-for="search-forward-input">
            This number (defaulted at 100) represents how many base pairs forward from the current start position 
            alternate start codons will be searched for.
            <b-form-input
              id="search-forward-input"
              type="number"
              v-model="forwardStartRange"
              required
              placeholder="Enter start search range"
            ></b-form-input>
          </b-form-group>
          <b-form-group label="Gap:" label-size="lg" label-for="gap-input">
            This number (defaulted at 10) represents the maximum acceptable number of base pairs between two adjacent genes. 
            Gaps greater than this number will be flagged.
            <b-form-input
              id="gap-input"
              type="number"
              v-model="gap"
              required
              placeholder="Enter acceptable gap length"
            ></b-form-input>
          </b-form-group>
          <b-form-group label="Overlap:" label-size="lg" label-for="overlap-input">
            This number (defaulted at 10) represents the maximum acceptable number of base pairs two adjacent genes overlap. 
            Overlaps greater than this number will be flagged.
            <b-form-input
              id="overlap-input"
              type="number"
              v-model="overlap"
              required
              placeholder="Enter acceptable overlap length."
            ></b-form-input>
          </b-form-group>
          <b-form-group label="Opposite Strand Gap:" label-size="lg" label-for="opposite-gap-input">
            This number (defaulted at 50) represents the minimum acceptable number of base pairs between two adjacent 
            genes on different strands. Gaps shorter than this number will be flagged.
            <b-form-input
              id="opposite-gap-input"
              type="number"
              v-model="oppositeGap"
              required
              placeholder="Enter acceptable gap length"
            ></b-form-input>
          </b-form-group>
          <b-form-group label="Minimum Gene Length:" label-size="lg" label-for="short-input">
            This number (defaulted at 200) represents the minimum acceptable number of base pairs in a gene. 
            Genes shorter than this number will be flagged.
            <b-form-input
              id="short-input"
              type="number"
              v-model="short"
              required
              placeholder="Enter acceptable minimum gene length"
            ></b-form-input>
          </b-form-group>
          <hr />
          <b-button type="submit" class="mt-3" block style="margin-top: 0px">
            <strong>Submit</strong>
          </b-button>
        </b-form>
      </b-modal>
  </div>
</template>

<script>
import axios from "axios";

export default {
  name: "Navbar",
  props: {
    upload: Boolean,
    // dnamaster: Boolean,
    blast: Boolean,
    annotations: Boolean,
    geneMap: Boolean,
    settings: Boolean,
    phageID: String,
  },

  data() {
    return {
      activeCLass: 'active',
      showSettings: false,
      gap: null,
      overlap: null,
      oppositeGap: null,
      backStartRange: null,
      forwardStartRange: null,
      short: null,
    }
  },

  methods: {
    getSettings() {
      this.showSettings = true,
      axios.get(process.env.VUE_APP_BASE_URL + `/settings/${this.$route.params.phageID}/none`)
      .then(response => {
        console.log(response.data)
        this.backStartRange = response.data.back_start_range;
        this.forwardStartRange = response.data.forward_start_range;
        this.gap = response.data.gap;
        this.overlap = response.data.overlap;
        this.oppositeGap = response.data.opposite_gap;
        this.short = response.data.short;
      })
    },

    onSubmit(evt) {
      evt.preventDefault();
      this.$refs.settingsModal.hide();
      var payload = this.gap + ',' + this.overlap + ',' + this.oppositeGap + ',' + this.backStartRange + ',' + this.forwardStartRange + ',' + this.short;
      axios.get(process.env.VUE_APP_BASE_URL + `/settings/${this.$route.params.phageID}/${payload}`)
        .then(response => {
          console.log(response.data);
        })
        .catch(error => {
          console.error(error);
        });
    },
  },
  
};
</script>

<style scoped>
.active {
  opacity: 1;
  visibility: visible;
  border-left-color: #4DB6AC;
  margin: 10px;
  transition: all 0.25s;
}
</style>