  <template>
    <div>
      <b-alert variant="success" show>
        {{ message }}
      </b-alert>
      <slot></slot>
      <br>
    </div>
  </template>

<script>
export default {
  props: ['message'],
};
</script>
