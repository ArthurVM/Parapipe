import os
from fpdf import FPDF

class PDF(FPDF):
    def __init__(self):
        super().__init__()
        self.WIDTH = 210
        self.HEIGHT = 297

    def write_subheader(self, text):
        """ write subheader text to pdf
        """
        self.set_font('Courier', 'B', 15)
        self.cell(self.WIDTH - 80)
        self.cell(-1, -1, text, 0, 0, 'R')
        self.ln(20)

    def write_title(self, sample_ID):

        self.set_font('Courier', 'B', 18)
        self.cell(self.WIDTH - 80)
        self.cell(-20, -10, f"Report for {sample_ID}", 0, 0, 'R')

    def header(self):
        # Custom logo and positioning
        # Create an `assets` folder and put any wide and short image inside
        # Name the image `logo.png`
        scriptdir = os.path.dirname(os.path.realpath(__file__))
        self.image(f'{scriptdir}/../resources/phw_logo.jpeg', 10, 8, 60)
        self.set_font('Courier', 'B', 11)
        self.cell(self.WIDTH - 80)
        self.cell(60, 1, 'Parapipe report', 0, 0, 'R')
        self.ln(20)

    def footer(self):
        # Page numbers in the footer
        self.set_y(-15)
        self.set_font('Courier', 'I', 8)
        self.set_text_color(128)
        self.cell(0, 10, 'Page ' + str(self.page_no()), 0, 0, 'C')

    def page_body(self, images, big_small=False, small_big=False):
        # Determine how many plots there are per page and set positions
        # and margins accordingly
        if big_small and small_big:
            raise ValueError("Key word arguments big_small and small_big may not be used together.")

        if len(images) == 3:
            self.image(images[0], 15, 25, self.WIDTH - 30)
            self.image(images[1], 15, self.WIDTH / 2 + 5, self.WIDTH - 30)
            self.image(images[2], 15, self.WIDTH / 2 + 90, self.WIDTH - 30)
        elif len(images) == 2:
            if big_small:
                self.image(images[0], 15, 25, self.WIDTH - 30)
                self.image(images[1], 15, self.WIDTH / 2 + 90, self.WIDTH - 30)
            elif small_big:
                self.image(images[0], 15, 25, self.WIDTH - 30)
                self.image(images[1], 15, self.WIDTH / 2 + 5, self.WIDTH - 30)
            else:
                self.image(images[0], 15, 25, self.WIDTH - 30)
                self.image(images[1], 15, self.WIDTH / 2 + 50, self.WIDTH - 30)
        else:
            self.image(images[0], 15, 40, self.WIDTH - 30)

    def print_page(self, images, big_small=False, small_big=False):
        # Generates the report
        self.add_page()
        self.page_body(images, big_small, small_big)
