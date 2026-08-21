use gpui::Context;
use settings::SettingsStore;
use std::time::{Duration, Instant};
use ui::App;

const SMOOTH_BLINK_FADE_DURATION: Duration = Duration::from_millis(150);

pub struct BlinkManager {
    blink_interval: Duration,
    blink_epoch: usize,
    /// Whether the blinking is paused.
    blinking_paused: bool,
    /// Whether the cursor should be visibly rendered or not.
    visible: bool,
    /// Whether the blinking is currently enabled.
    enabled: bool,
    /// Whether the blinking is enabled in the settings.
    blink_enabled_in_settings: fn(&App) -> bool,
    /// Whether smooth blink transitions are enabled.
    smooth_blink_enabled: bool,
    /// Opacity the current fade started from.
    fade_from: f32,
    /// Opacity the current fade ends at.
    fade_target: f32,
    /// When the current fade started; `None` when settled at `fade_target`.
    fade_start: Option<Instant>,
}

/// Ease a fade progress value in `[0.0, 1.0]` with a smoothstep curve.
fn ease_smoothstep(progress: f32) -> f32 {
    let progress = progress.clamp(0.0, 1.0);
    progress * progress * (3.0 - 2.0 * progress)
}

impl BlinkManager {
    pub fn new(
        blink_interval: Duration,
        blink_enabled_in_settings: fn(&App) -> bool,
        cx: &mut Context<Self>,
    ) -> Self {
        // Make sure we blink the cursors if the setting is re-enabled
        cx.observe_global::<SettingsStore>(move |this, cx| {
            this.blink_cursors(this.blink_epoch, cx)
        })
        .detach();

        Self {
            blink_interval,
            blink_epoch: 0,
            blinking_paused: false,
            visible: true,
            enabled: false,
            blink_enabled_in_settings,
            smooth_blink_enabled: false,
            fade_from: 1.0,
            fade_target: 1.0,
            fade_start: None,
        }
    }

    pub fn set_smooth_blink_enabled(&mut self, enabled: bool) {
        self.smooth_blink_enabled = enabled;
        self.settle_opacity(if self.visible { 1.0 } else { 0.0 });
    }

    /// Jump to the given opacity with no fade.
    fn settle_opacity(&mut self, opacity: f32) {
        self.fade_from = opacity;
        self.fade_target = opacity;
        self.fade_start = None;
    }

    /// Start a finite fade from the currently displayed opacity to `target`.
    ///
    /// With smooth blink disabled (or when already at `target`) this settles
    /// instantly, so no animation frames are requested.
    fn start_fade(&mut self, target: f32) {
        if !self.smooth_blink_enabled {
            self.settle_opacity(target);
            return;
        }
        let current = self.smooth_opacity();
        if (current - target).abs() < f32::EPSILON {
            self.settle_opacity(target);
            return;
        }
        self.fade_from = current;
        self.fade_target = target;
        self.fade_start = Some(Instant::now());
    }

    /// The smooth-blink opacity at this instant, derived purely from the fade
    /// start time. The fade is a finite ramp: it reaches `fade_target` exactly
    /// when `SMOOTH_BLINK_FADE_DURATION` elapses and never animates past it.
    fn smooth_opacity(&self) -> f32 {
        match self.fade_start {
            None => self.fade_target,
            Some(start) => {
                let progress = start.elapsed().as_secs_f32()
                    / SMOOTH_BLINK_FADE_DURATION.as_secs_f32();
                self.fade_from + (self.fade_target - self.fade_from) * ease_smoothstep(progress)
            }
        }
    }

    /// The smooth-blink opacity (0.0 = hidden, 1.0 = solid).
    ///
    /// This deliberately does not call `cx.notify()`: redraws during a fade are
    /// driven by the cursor animation-frame loop, and notifying here would
    /// invalidate the window's cached scene every frame.
    pub fn opacity(&self) -> f32 {
        if !self.smooth_blink_enabled {
            return if self.visible { 1.0 } else { 0.0 };
        }
        self.smooth_opacity()
    }

    pub fn is_smooth_blink_animating(&self) -> bool {
        self.smooth_blink_enabled
            && self
                .fade_start
                .is_some_and(|start| start.elapsed() < SMOOTH_BLINK_FADE_DURATION)
    }

    fn next_blink_epoch(&mut self) -> usize {
        self.blink_epoch += 1;
        self.blink_epoch
    }

    /// Show the cursor immediately and reschedule blinking. Called when the
    /// user moves the cursor or types; blinking resumes after a short delay.
    pub fn pause_blinking(&mut self, cx: &mut Context<Self>) {
        self.show_cursor(cx);
        if self.smooth_blink_enabled {
            self.blinking_paused = true;
        }

        let epoch = self.next_blink_epoch();
        let interval = Duration::from_millis(500);
        cx.spawn(async move |this, cx| {
            cx.background_executor().timer(interval).await;
            this.update(cx, |this, cx| this.resume_cursor_blinking(epoch, cx))
        })
        .detach();
    }

    fn resume_cursor_blinking(&mut self, epoch: usize, cx: &mut Context<Self>) {
        if epoch == self.blink_epoch {
            self.blinking_paused = false;
            if self.smooth_blink_enabled {
                // Show the cursor solid for one interval before the next
                // toggle, so resuming does not fade out immediately.
                self.show_cursor(cx);
                let interval = self.blink_interval;
                cx.spawn(async move |this, cx| {
                    cx.background_executor().timer(interval).await;
                    if let Some(this) = this.upgrade() {
                        this.update(cx, |this, cx| this.blink_cursors(epoch, cx));
                    }
                })
                .detach();
            } else {
                self.blink_cursors(epoch, cx);
            }
        }
    }

    fn blink_cursors(&mut self, epoch: usize, cx: &mut Context<Self>) {
        if (self.blink_enabled_in_settings)(cx) {
            if epoch == self.blink_epoch && self.enabled && !self.blinking_paused {
                self.visible = !self.visible;
                self.start_fade(if self.visible { 1.0 } else { 0.0 });
                cx.notify();

                let epoch = self.next_blink_epoch();
                let interval = self.blink_interval;
                cx.spawn(async move |this, cx| {
                    cx.background_executor().timer(interval).await;
                    if let Some(this) = this.upgrade() {
                        this.update(cx, |this, cx| this.blink_cursors(epoch, cx));
                    }
                })
                .detach();
            }
        } else {
            self.show_cursor(cx);
        }
    }

    pub fn show_cursor(&mut self, cx: &mut Context<BlinkManager>) {
        if self.smooth_blink_enabled {
            self.visible = true;
            self.settle_opacity(1.0);
            cx.notify();
        } else if !self.visible {
            self.visible = true;
            cx.notify();
        }
    }

    /// Enable the blinking of the cursor.
    pub fn enable(&mut self, cx: &mut Context<Self>) {
        if self.enabled {
            return;
        }

        self.enabled = true;
        // Set cursors as invisible and start blinking: this causes cursors
        // to be visible during the next render.
        self.visible = false;
        self.start_fade(0.0);
        self.blink_cursors(self.blink_epoch, cx);
    }

    /// Disable the blinking of the cursor.
    pub fn disable(&mut self, _cx: &mut Context<Self>) {
        self.visible = false;
        self.enabled = false;
        self.start_fade(0.0);
    }

    pub fn visible(&self) -> bool {
        self.visible
    }

    #[cfg(test)]
    pub(crate) fn enabled(&self) -> bool {
        self.enabled
    }

    pub fn should_render(&self) -> bool {
        if self.smooth_blink_enabled {
            self.visible || self.is_smooth_blink_animating()
        } else {
            self.visible
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gpui::{AppContext, TestAppContext};

    #[gpui::test]
    async fn cursor_move_forces_cursor_visible_immediately(cx: &mut TestAppContext) {
        let blink_manager =
            cx.new(|cx| BlinkManager::new(Duration::from_millis(500), |_| true, cx));

        blink_manager.update(cx, |blink_manager: &mut BlinkManager, cx| {
            blink_manager.disable(cx);
            blink_manager.set_smooth_blink_enabled(true);

            assert_eq!(blink_manager.opacity(), 0.0);
            blink_manager.pause_blinking(cx);
            assert_eq!(blink_manager.opacity(), 1.0);
            assert!(blink_manager.should_render());
            assert!(!blink_manager.is_smooth_blink_animating());
        });
    }

    #[gpui::test]
    async fn smooth_blink_fade_is_a_finite_ramp(cx: &mut TestAppContext) {
        let blink_manager =
            cx.new(|cx| BlinkManager::new(Duration::from_millis(500), |_| true, cx));

        blink_manager.update(cx, |blink_manager: &mut BlinkManager, _cx| {
            blink_manager.set_smooth_blink_enabled(true);

            assert_eq!(blink_manager.opacity(), 1.0);
            blink_manager.start_fade(0.0);
            assert!(blink_manager.is_smooth_blink_animating());
        });

        // Real wall-clock wait: the fade is driven by `Instant`, not the test
        // executor. After the 150ms ramp the fade must be exactly settled so
        // animation-only frames stop.
        std::thread::sleep(SMOOTH_BLINK_FADE_DURATION + Duration::from_millis(50));

        blink_manager.update(cx, |blink_manager: &mut BlinkManager, _cx| {
            assert!(!blink_manager.is_smooth_blink_animating());
            assert_eq!(blink_manager.opacity(), 0.0);
        });
    }

    #[test]
    fn smoothstep_ease_hits_endpoints_exactly() {
        assert_eq!(ease_smoothstep(0.0), 0.0);
        assert_eq!(ease_smoothstep(1.0), 1.0);
        // Values past the ramp clamp to the endpoint instead of overshooting.
        assert_eq!(ease_smoothstep(2.5), 1.0);
        assert_eq!(ease_smoothstep(-1.0), 0.0);
        // Monotonic in between.
        assert!(ease_smoothstep(0.25) < ease_smoothstep(0.5));
        assert!(ease_smoothstep(0.5) < ease_smoothstep(0.75));
    }
}
