import logging
from joblib import delayed, Parallel
from obspy import UTCDateTime
from obspy.signal.trigger import ar_pick, pk_baer
from obspy.core.event import Event, Pick, WaveformStreamID, Amplitude
from obspy.core.event import CreationInfo, Comment, Origin

from phasepapy.phasepicker.aicdpicker import AICDPicker
from phasepapy.phasepicker.ktpicker import KTPicker
from phasepapy.phasepicker.fbpicker import FBPicker

PST_AUTHOR = ':GA-PST'
AGENCY_URI = 'GA'
log = logging.getLogger(__name__)


def _naive_phase_type(tr):
    """
    copied from EQcorrscan.git
    :param tr: obspy.trace object
    :return: str
    """
    # We are going to assume, for now, that if the pick is made on the
    # horizontal channel then it is an S, otherwise we will assume it is
    # a P-phase: obviously a bad assumption...
    if tr.stats.channel[-1] == 'Z':
        phase = 'P'
    else:
        phase = 'S'
    return phase


# write custom picker classes here
class PKBaer:
    """
    our wrapper class for obspy pk_baer picker class.

    For inputs and outputs refer to obspy.signal.trigge.pk_bear class
    """

    def __init__(self, tdownmax=20, tupevent=60, thr1=7.0, thr2=12.0,
                 preset_len=100, p_dur=100):
        self.tdownmax = tdownmax
        self.tupevent = tupevent
        self.thr1 = thr1
        self.thr2 = thr2
        self.preset_len = preset_len
        self.p_dur = p_dur

    def picks(self, tr):
        """
        Parameters
        ----------
        tr: obspy.core.trace.Trace
            obspy trace instance
        samp_int: int
            number of samples per second

        Returns
        -------
        pptime: int
            pptime sample number of parrival
        pptime sample number of parrival
        pfm: str
            pfm direction of first motion (U or D)
        """
        reltrc = tr.data
        samp_int = tr.stats.sampling_rate
        pptime, pfm = pk_baer(reltrc=reltrc,
                              samp_int=samp_int,
                              tdownmax=self.tdownmax,
                              tupevent=self.tupevent,
                              thr1=self.thr1,
                              thr2=self.thr2,
                              preset_len=self.preset_len,
                              p_dur=self.p_dur)

        return pptime / samp_int, pfm


class AICDPickerGA(AICDPicker):
    """
    Our adaptation of AICDpicker.
    See docs for AICDpicker.
    """

    def __init__(self, t_ma=3, nsigma=6, t_up=0.2, nr_len=2,
                 nr_coeff=2, pol_len=10, pol_coeff=10,
                 uncert_coeff=3):
        AICDPicker.__init__(self,
                            t_ma=t_ma,
                            nsigma=nsigma,
                            t_up=t_up,
                            nr_len=nr_len,
                            nr_coeff=nr_coeff,
                            pol_len=pol_len,
                            pol_coeff=pol_coeff,
                            uncert_coeff=uncert_coeff)

    def event(self, st, config):
        event = Event()
        # event.origins.append(Origin())
        creation_info = CreationInfo(author=self.__class__.__name__ +
                                            PST_AUTHOR,
                                     creation_time=UTCDateTime(),
                                     agency_uri=AGENCY_URI)

        event.creation_info = creation_info
        event.comments.append(Comment(text='pst-aicdpicker'))

        # note Parallel might cause issues during mpi processing, ok for now
        # TODO: should split this to multiple MPI processes
        res = Parallel(n_jobs=-1)(delayed(self._pick_parallel)(tr, config)
                                  for tr in st)

        for pks, ws, ps in res:
            for p in pks:
                event.picks.append(
                    Pick(waveform_id=ws, phase_hint=ps, time=p,
                         creation_info=creation_info)
                    # FIXME: same creation info for all picks
                )
        return event

    def _pick_parallel(self, tr, config):
        if config.detrend:
            tr = tr.detrend(config.detrend)
            log.info("Applied '{}' detrend to trace".format(config.detrend))
        if config.filter:
            tr = tr.filter(config.filter['type'], ** config.filter['params'])
            log.info("Applied '{}' filter with params: {}".format(
                config.filter['type'], config.filter['params']))
        _, picks, polarity, snr, uncertainty = self.picks(tr=tr)
        phase = _naive_phase_type(tr)
        wav_id = WaveformStreamID(station_code=tr.stats.station,
                                  channel_code=tr.stats.channel,
                                  network_code=tr.stats.network,
                                  location_code=tr.stats.location
                                  )
        return picks, wav_id, phase


pickermaps = {
    'aicdpicker': AICDPickerGA,
    'fbpicker': FBPicker,
    'ktpicker': KTPicker,
    # 'pkbaer': PKBaer,
}
