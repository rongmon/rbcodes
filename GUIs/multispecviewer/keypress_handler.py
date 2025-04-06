from PyQt5.QtCore import Qt

class KeyPressHandler:
    def __init__(self, plot_canvas):
        """Initializes the keypress handler with a Matplotlib canvas."""
        self.canvas = plot_canvas.canvas
        self.figure = plot_canvas.figure
        self.ax_list = self.figure.get_axes()

        # Connect Matplotlib event (for mouse focus)
        self.canvas.mpl_connect("key_press_event", self.on_key_press)

    def on_key_press(self, event):
        """Handles both Matplotlib and Qt keypress events."""
        key = event.key if hasattr(event, "key") else event.text()
        print(f"Key pressed: {key}")  # Debugging print

        if event.inaxes:  # Check if event is inside an axis
            print(f"Event Coordinates - xdata: {event.xdata}, ydata: {event.ydata}")
        else:
            print("Keypress ignored: Not inside an axis.")  # Event occurred outside axis

        for ax in self.ax_list:
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()

            if key == 'x' and event.xdata is not None:
                print(f"Setting min x-limit to {event.xdata}")
                ax.set_xlim([event.xdata, xlim[1]])
            elif key == 'X' and event.xdata is not None:
                print(f"Setting max x-limit to {event.xdata}")
                ax.set_xlim([xlim[0], event.xdata])
            elif key == 't' and event.ydata is not None:
                print(f"Setting max y-limit to {event.ydata}")
                ax.set_ylim([ylim[0], event.ydata])
            elif key == 'b' and event.ydata is not None:
                print(f"Setting min y-limit to {event.ydata}")
                ax.set_ylim([event.ydata, ylim[1]])

        self.canvas.draw_idle()  # Ensure UI updates immediately
